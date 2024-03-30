#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <iostream>
#include <string.h>
#include <math.h>
#include <chrono>
#include <omp.h>

int GuidedBilateralFilterStep(int dimx, int dimy, int ncol, unsigned char const *orig, unsigned char const *guide, int demisize,
							  float sscale, float iscale, float ipower, float gscale, float gpower, float *filtered)
{
	float *sweight = NULL, iweight[257], gweight[256];

	/* spatial weight */
	if ((sweight = (float *)malloc((demisize + 1) * sizeof(float))) == NULL)
		return (0);
	for (int i = 0; i <= demisize; i++)
	{
		if (sscale > 0.0f)
			sweight[i] = exp(-0.5f * (float)(i * i) / (sscale * sscale));
		else
			sweight[i] = 1.0f;
	}

	/* intensity weight */
	for (int i = 0; i <= 256; i++)
	{
		if (ipower != 1.0f)
			iweight[i] = pow(1.0f + (float)(i * i) / (iscale * iscale), ipower - 1.0f);
		else
			iweight[i] = 1.0f;
	}

	/* guide weight */
	for (int i = 0; i <= 255; i++)
	{
		if (gpower != 0.0f)
			gweight[i] = exp(-(pow(1.0f + (float)(i * i) / (gscale * gscale), gpower) - 1.0f) / gpower);
		else
			gweight[i] = 1.0f / (1.0f + (float)(i * i) / (gscale * gscale));
	}

	// this loop is the slow part
	#pragma omp parallel for collapse(2)
	for (int j = 0; j < dimy; j++)
	{
		for (int i = 0; i < dimx; i++)
		{
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
	}

	free(sweight);

	return (1);
}

int GuidedBilateralFilter(int dimx, int dimy, int ncol, unsigned char const *orig, unsigned char const *guide, int demisize, float sscale, float iscale, float ipower, float gscale, float gpower, unsigned char *result)
{
	int i, num = 8;

	/* alloc */
	float *filtered = new float[dimx * dimy];

	/* init image */
	for (i = 0; i < dimx * dimy; i++)
		filtered[i] = (float)(orig[i]);

	/* GNC */
	if (ipower <= 1.0f)
	{
		if (!GuidedBilateralFilterStep(dimx, dimy, ncol, orig, guide, demisize, 0.0, iscale, 1.0, gscale * 5.0, gpower, filtered))
			return (0);
		num--;
	}

	if (ipower <= 0.5f)
	{
		if (!GuidedBilateralFilterStep(dimx, dimy, ncol, orig, guide, demisize, sscale, iscale, 0.5, gscale, gpower, filtered))
			return (0);
		num--;
	}

	if (ipower <= 0.0f)
	{
		if (!GuidedBilateralFilterStep(dimx, dimy, ncol, orig, guide, demisize, sscale, iscale, 0.0, gscale, gpower, filtered))
			return (0);
		num--;
	}

	/* final */
	for (i = 0; i < num; i++)
	{
		if (!GuidedBilateralFilterStep(dimx, dimy, ncol, orig, guide, demisize, sscale, iscale, ipower, gscale, gpower, filtered))
			return (0);
	}

	for (i = 0; i < dimx * dimy; i++)
		result[i] = (unsigned char)(filtered[i]);

	delete[] filtered;

	return (1);
}

// very slow implementation, to improve
// - use cv::cuda functions
// - do not split and merge the color channels, change the above functions for the images with stacked color channels
// - filter is parallelizable, implement it as a cuda plugin
// - decrease the iteration count num = 8 in GuidedBilateralFilter()

// i tried to add cpu multithreading. 647ms in karagag server.
cv::Mat GuidedBilateralFilterToCVImage(cv::Mat origimg_, cv::Mat guideimg_){
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


	omp_set_num_threads(32);

	// cv::imshow("orig", origimg_);
	// cv::imshow("guide", guideimg_);

	cv::Mat origimg[3], guideimg[3];
	cv::split(origimg_, origimg);
	cv::split(guideimg_, guideimg);

	std::vector<cv::Mat> resultmatIIminusIJ;
	resultmatIIminusIJ.reserve(3);

	// bgr color channels loop
	// TODO: i tried to parallize here, but could not
	for (int i = 0; i < 3; i++)
	{
		unsigned char * resultIJ = new unsigned char[origimg[i].rows * origimg[i].cols];
		unsigned char * resultII = new unsigned char[origimg[i].rows * origimg[i].cols];

		GuidedBilateralFilter(origimg[i].rows, origimg[i].cols, origimg[i].channels(), origimg[i].data, guideimg[i].data, hwsize, sscale, iscale, ipower, gscale, gpower, resultIJ);
		cv::Mat resultmatIJ = cv::Mat(origimg[i].rows, origimg[i].cols, CV_8U, resultIJ);

		GuidedBilateralFilter(origimg[i].rows, origimg[i].cols, origimg[i].channels(), origimg[i].data, origimg[i].data, hwsize, sscale, iscale, ipower, gscale, gpower, resultII);
		cv::Mat resultmatII = cv::Mat(origimg[i].rows, origimg[i].cols, CV_8U, resultII);

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

int main()
{
	cv::Mat origimg_ = cv::imread("../input_images/makale_1.png", cv::IMREAD_COLOR);
	origimg_.convertTo(origimg_, CV_8U); // just for safety
	cv::Mat guideimg_ = cv::imread("../input_images/makale_0.png", cv::IMREAD_COLOR);
	guideimg_.convertTo(guideimg_, CV_8U); // just for safety

	int n_iter = 1;
	auto start = std::chrono::steady_clock::now();
	cv::Mat result;
	for(int i = 0; i < n_iter; i++)
	result = GuidedBilateralFilterToCVImage(origimg_, guideimg_);    
    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time in milliseconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/n_iter
        << " ms\n";

	// cv::imshow("distance all channels", result);
	// int k = cv::waitKey(0);
	cv::imwrite("../output_images/result_cpu.png", result);

	return 0;
}