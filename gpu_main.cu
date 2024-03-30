#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <iostream>
#include <string.h>
#include <cmath>
#include <chrono>
#include <map>

__global__ void bilateralKernel(int dimx, int dimy, int ncol, unsigned char *orig, unsigned char *guide, int demisize,
								float *sweight, float *iweight, float *gweight,
								float *filtered)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;

	if (j >= dimy || i >= dimx)
		return;

	int value, ediff, currentGuide[3], diffGuide;
	float wguide, somme, poids, pixelMoy, currentIntensity, diff, rdiff;
	somme = 1e-6f;
	pixelMoy = 0.0f;
	currentIntensity = filtered[j * dimx + i];
	currentGuide[0] = guide[j * dimx + i];
	if (ncol == 3)
	{
		currentGuide[1] = guide[dimx * dimy + j * dimx + i];
		currentGuide[2] = guide[2 * dimx * dimy + j * dimx + i];
	}
	// don't need to parallize here since only 2x2 max
	for (int k = -demisize; k <= demisize; k++)
	{
		if ((j + k >= 0) && (j + k < dimy))
		{
			for (int l = -demisize; l <= demisize; l++)
			{
				if ((i + l >= 0) && (i + l < dimx))
				{
					value = orig[(j + k) * dimx + i + l];
					diff = fabs((float)value - currentIntensity);
					ediff = (int)floor(diff);
					rdiff = diff - (float)ediff;
					diffGuide = abs(guide[(j + k) * dimx + i + l] - currentGuide[0]);
					wguide = gweight[diffGuide];
					if (ncol == 3)
					{
						diffGuide = abs(guide[dimx * dimy + (j + k) * dimx + i + l] - currentGuide[1]);
						wguide *= gweight[diffGuide];
						diffGuide = abs(guide[2 * dimx * dimy + (j + k) * dimx + i + l] - currentGuide[2]);
						wguide *= gweight[diffGuide];
					}
					poids = ((1.0f - rdiff) * iweight[ediff] + rdiff * iweight[ediff + 1]) * sweight[abs(k)] * sweight[abs(l)] * wguide;
					somme += poids;
					pixelMoy += poids * (float)value;
				}
			}
		}
	}

	filtered[j * dimx + i] = pixelMoy / somme;
}

class GuidedBilateralFilterGPU
{
public:
	// Guided Bilateral Filter parameters
	int hwsize = 2;
	float sscale = 1.5f, iscale = 10.0f, ipower = 0.0f, gscale = 10.0f, gpower = 1.0f;

	// Threshold parameter
	int threshold = 80;

	// Opening parameters
	int morph_size = 1;
	cv::Mat element = getStructuringElement(
		cv::MORPH_ELLIPSE,
		cv::Size(2 * morph_size + 1,
				 2 * morph_size + 1),
		cv::Point(morph_size,
				  morph_size));

	cv::Mat origimg[3], guideimg[3];

	float *filtered_d;
	unsigned char *orig_d;
	unsigned char *guide_d;

	float *sweight, *iweight, *gweight;
	float *sweight_d, *iweight_d, *gweight_d;

	float *filtered_cpu;
	int size_, size;

	GuidedBilateralFilterGPU(int rows, int cols)
	{
		size_ = rows * cols;
		size = size_ * sizeof(float);

		cudaMalloc((float **)&filtered_d, size);
		cudaMalloc((unsigned char **)&orig_d, size_);
		cudaMalloc((unsigned char **)&guide_d, size_);

		sweight = (float *)malloc((hwsize + 1) * sizeof(float));
		cudaMalloc((float **)&sweight_d, (hwsize + 1) * sizeof(float));
		cudaMalloc((float **)&iweight_d, 257 * sizeof(float));
		cudaMalloc((float **)&gweight_d, 256 * sizeof(float));

		filtered_cpu = (float *)malloc(size);
	}

	std::map<std::pair<float, float>, float *> iweights;
	void iweightcalculation(float iscale, float ipower)
	{
		auto ii = iweights.find(std::make_pair(iscale, ipower));
		if (ii != iweights.end())
		{
			iweight = ii->second;
		}
		else
		{
			float *new_iweight = (float *)malloc(257 * sizeof(float));
			/* intensity weight */
			for (int i = 0; i <= 256; i++)
			{
				if (ipower != 1.0f)
					new_iweight[i] = pow(1.0f + (float)(i * i) / (iscale * iscale), ipower - 1.0f);
				else
					new_iweight[i] = 1.0f;
			}
			iweights[std::make_pair(iscale, ipower)] = new_iweight;
			iweight = new_iweight;
		}
	}

	std::map<std::pair<float, float>, float *> gweights;
	void gweightcalculation(float gscale, float gpower)
	{
		auto ii = gweights.find(std::make_pair(gscale, gpower));
		if (ii != gweights.end())
		{
			gweight = ii->second;
		}
		else
		{
			float *new_gweight = (float *)malloc(256 * sizeof(float));
			/* guide weight */
			for (int i = 0; i <= 255; i++)
			{
				if (gpower != 0.0f)
					new_gweight[i] = exp(-(pow(1.0f + (float)(i * i) / (gscale * gscale), gpower) - 1.0f) / gpower);
				else
					new_gweight[i] = 1.0f / (1.0f + (float)(i * i) / (gscale * gscale));
			}
			gweights[std::make_pair(gscale, gpower)] = new_gweight;
			gweight = new_gweight;
		}
	}

	int GuidedBilateralFilterStep(int dimx, int dimy, int ncol, unsigned char *orig, unsigned char *guide, int demisize,
								  float sscale, float iscale, float ipower, float gscale, float gpower)
	{
		for (int i = 0; i <= demisize; i++)
		{
			if (sscale > 0.0f)
				sweight[i] = exp(-0.5f * (float)(i * i) / (sscale * sscale));
			else
				sweight[i] = 1.0f;
		}

		iweightcalculation(iscale, ipower);
		gweightcalculation(gscale, gpower);

		cudaMemcpy(sweight_d, sweight, (demisize + 1) * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(iweight_d, iweight, 257 * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(gweight_d, gweight, 256 * sizeof(float), cudaMemcpyHostToDevice);

		dim3 block(16, 16);
		dim3 grid((dimx + 15) / 16, (dimy + 15) / 16);
		bilateralKernel<<<grid, block>>>(dimx, dimy, ncol, orig, guide, demisize,
										 sweight_d, iweight_d, gweight_d,
										 filtered_d);

		return (1);
	}

	int GuidedBilateralFilter(int dimx, int dimy, int ncol, unsigned char *orig, unsigned char *guide, int demisize, float sscale, float iscale, float ipower, float gscale, float gpower, unsigned char *result)
	{
		int i, num = 8;

		/* init image */
		for (i = 0; i < dimx * dimy; i++)
			filtered_cpu[i] = (float)(orig[i]);

		cudaMemcpy(filtered_d, filtered_cpu, size, cudaMemcpyHostToDevice);
		cudaMemcpy(orig_d, orig, (dimx * dimy), cudaMemcpyHostToDevice);
		cudaMemcpy(guide_d, guide, (dimx * dimy), cudaMemcpyHostToDevice);

		/* GNC */
		if (ipower <= 1.0f)
		{
			if (!GuidedBilateralFilterStep(dimx, dimy, ncol, orig_d, guide_d, demisize, 0.0, iscale, 1.0, gscale * 5.0, gpower))
				return (0);
			num--;
		}

		if (ipower <= 0.5f)
		{
			if (!GuidedBilateralFilterStep(dimx, dimy, ncol, orig_d, guide_d, demisize, sscale, iscale, 0.5, gscale, gpower))
				return (0);
			num--;
		}

		if (ipower <= 0.0f)
		{
			if (!GuidedBilateralFilterStep(dimx, dimy, ncol, orig_d, guide_d, demisize, sscale, iscale, 0.0, gscale, gpower))
				return (0);
			num--;
		}

		/* final */
		for (i = 0; i < num; i++)
		{
			if (!GuidedBilateralFilterStep(dimx, dimy, ncol, orig_d, guide_d, demisize, sscale, iscale, ipower, gscale, gpower))
				return (0);
		}

		cudaMemcpy(filtered_cpu, filtered_d, size, cudaMemcpyDeviceToHost);

		for (i = 0; i < dimx * dimy; i++)
			result[i] = (unsigned char)(filtered_cpu[i]);

		// cudaError_t error_check = cudaGetLastError();printf("%s\n", cudaGetErrorString(error_check));

		return (1);
	}

	// very slow implementation, to improve
	// - use cv::cuda functions
	// - do not split and merge the color channels, change the above functions for the images with stacked color channels
	cv::Mat Execute(cv::Mat origimg_, cv::Mat guideimg_)
	{
		// cv::imshow("orig", origimg_);
		// cv::imshow("guide", guideimg_);

		cv::split(origimg_, origimg);
		cv::split(guideimg_, guideimg);

		std::vector<cv::Mat> resultmatIIminusIJ;
		resultmatIIminusIJ.reserve(3);

		// bgr color channels loop
		// TODO: i tried to parallize here, but could not
		for (int i = 0; i < 3; i++)
		{
			unsigned char *resultIJ = new unsigned char[origimg[i].rows * origimg[i].cols];
			unsigned char *resultII = new unsigned char[origimg[i].rows * origimg[i].cols];

			GuidedBilateralFilter(origimg[i].rows, origimg[i].cols, origimg[i].channels(), origimg[i].data, guideimg[i].data, hwsize, sscale, iscale, ipower, gscale, gpower, resultIJ);
			cv::Mat resultmatIJ = cv::Mat(origimg[i].rows, origimg[i].cols, CV_8U, resultIJ);

			GuidedBilateralFilter(origimg[i].rows, origimg[i].cols, origimg[i].channels(), origimg[i].data, origimg[i].data, hwsize, sscale, iscale, ipower, gscale, gpower, resultII);
			cv::Mat resultmatII = cv::Mat(origimg[i].rows, origimg[i].cols, CV_8U, resultII);
			// cv::imwrite("../output_images/result_gpu_IJ.png", resultmatIJ);
			// cv::imwrite("../output_images/result_gpu_II.png", resultmatII);

			cv::Mat resultmatIIminusIJ_channel;
			cv::absdiff(resultmatII, resultmatIJ, resultmatIIminusIJ_channel);

			cv::threshold(resultmatIIminusIJ_channel, resultmatIIminusIJ_channel, threshold, 255, 1);

			morphologyEx(resultmatIIminusIJ_channel, resultmatIIminusIJ_channel,
						 cv::MORPH_OPEN, element,
						 cv::Point(-1, -1), 2);

			resultmatIIminusIJ.emplace_back(resultmatIIminusIJ_channel);

			// cv::imshow("resIJ channel:" + std::to_string(i), resultmatIJ);
			// cv::imshow("resII channel:" + std::to_string(i), resultmatII);
			// cv::imshow("distance channel:" + std::to_string(i), resultmatIIminusIJ_channel);
		}

		cv::Mat mergedresultmatIIminusIJ;
		merge(resultmatIIminusIJ, mergedresultmatIIminusIJ);

		return mergedresultmatIIminusIJ;
	}

	~GuidedBilateralFilterGPU()
	{
		cudaFree(filtered_d);
		cudaFree(orig_d);
		cudaFree(guide_d);

		cudaFree(sweight_d);
		cudaFree(iweight_d);
		cudaFree(gweight_d);

		free(filtered_cpu);

		for(auto ii : iweights) free(ii.second);
		for(auto ii : gweights) free(ii.second);
		free(sweight);
	}
};

int main()
{
	cv::Mat origimg_ = cv::imread("../input_images/makale_1.png", cv::IMREAD_COLOR);
	origimg_.convertTo(origimg_, CV_8U); // just for safety
	cv::Mat guideimg_ = cv::imread("../input_images/makale_0.png", cv::IMREAD_COLOR);
	guideimg_.convertTo(guideimg_, CV_8U); // just for safety

	GuidedBilateralFilterGPU gbFilter(origimg_.rows, origimg_.cols);

	int n_iter = 1;
	auto start = std::chrono::steady_clock::now();
	cv::Mat result;
	for (int i = 0; i < n_iter; i++)
		result = gbFilter.Execute(origimg_, guideimg_);
	auto end = std::chrono::steady_clock::now();
	std::cout << "Elapsed time in milliseconds for one frame: "
			  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / n_iter
			  << " ms\n";

	// cv::imshow("distance all channels", result);
	// int k = cv::waitKey(0);
	cv::imwrite("../output_images/result_gpu.png", result);

	return 0;
}