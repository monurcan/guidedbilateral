/*
 * =============================================================
 * jointcrossbilateralfilter_mex.c 
 *
 * Matlab-mex interface to the joint/cross bilateral filtering 
 * of an image
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2015 IFSTTAR
 * 
 * Jean-Philippe TAREL (16 february 2015)
 *
 * =============================================================
 */

/* $Revision: 0.0.0.1 $ */

#include <string.h>
#include <math.h>

#include "mex.h"

int JointCrossBilateralFilter(int dimx, int dimy, int ncol, unsigned char *orig, unsigned char *guide, int demisize, 
	float sscale, float gscale, unsigned char *result)
{
 int i, j, k, l, value, currentGuide[3], diffGuide;
 float somme, poids, pixelMoy, wguide, *sweight=NULL, gweight[256];    

 /* spatial weight */
 if ( (sweight = (float *)malloc( (demisize+1) * sizeof(float) ) ) == NULL) return(0);
 for(i=0; i<=demisize; i++) {
	if (sscale>0.0f) sweight[i]=exp(-0.5f*((float)i*(float)i)/(sscale*sscale));
	else sweight[i]=1.0f;
 }

 /* guide weight */
 for(i=0; i<=255; i++) gweight[i]=exp(-0.5f*(float)(i*i)/(gscale*gscale));

 for(j=0; j<dimy; j++) {
 	for(i=0; i<dimx; i++) {
		somme = 1e-6f;
		pixelMoy = 0.0f;
		currentGuide[0]=guide[j*dimx+i];
		if (ncol==3) {
			currentGuide[1]=guide[dimx*dimy+j*dimx+i];
			currentGuide[2]=guide[2*dimx*dimy+j*dimx+i];
		}
		for(k=-demisize; k<=demisize; k++) {  
			if ((j+k>=0) && (j+k<dimy)) {
				for(l=-demisize; l<=demisize; l++) {
			  		if ((i+l>=0) && (i+l<dimx)) {
			    			value = orig[(j+k)*dimx+i+l];
			    			diffGuide = abs(guide[(j+k)*dimx+i+l]-currentGuide[0]);
						wguide=gweight[diffGuide];
						if (ncol==3) {
							diffGuide = abs(guide[dimx*dimy+(j+k)*dimx+i+l]-currentGuide[1]);
							wguide*=gweight[diffGuide];
							diffGuide = abs(guide[2*dimx*dimy+(j+k)*dimx+i+l]-currentGuide[2]);
							wguide*=gweight[diffGuide];
						}
						poids = sweight[abs(k)]*sweight[abs(l)]*wguide;
			    			somme += poids;
			    			pixelMoy += poids*(float)value;
					}
			  	}
			}
		}
		poids = pixelMoy/somme;
		result[j*dimx+i]=(unsigned char)(poids);
	}
 }
 free(sweight);

 return(1);
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  unsigned char *input=NULL, *guide=NULL, *result=NULL;
  int ndims, ndimsGuide, ncol, dims[3], hwsize=2;
  float sscale=1.5f, gscale=10.0f;

  /*  Check for proper number of arguments. */
  if (nrhs < 2) mexErrMsgTxt("Two input argument are at least required.");

  /* Check to make sure the first input argument is an uint8 image */
  if (!mxIsClass(prhs[0],"uint8")) mexErrMsgTxt("Input image first should be uint8.");
  if (!mxIsClass(prhs[1],"uint8")) mexErrMsgTxt("Input image second should be uint8.");

  /* Create a pointer to the input image */
  input = (unsigned char *)mxGetData(prhs[0]);
  guide = (unsigned char *)mxGetData(prhs[1]);

  if (nrhs > 2) hwsize=(int)mxGetScalar(prhs[2]);
  if (nrhs > 3) sscale=(float)mxGetScalar(prhs[3]);
  if (nrhs > 4) gscale=(float)mxGetScalar(prhs[4]);
  if (nrhs > 5) mexErrMsgTxt("Too many input arguments.");
  if (nlhs < 1) mexErrMsgTxt("One output argument required.");
  if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");

  /* Get the dimensions of the image input */
  ndims = mxGetNumberOfDimensions(prhs[0]);
  /* hauteur image */
  dims[0] = mxGetDimensions(prhs[0])[0];
  if (dims[0] != mxGetDimensions(prhs[1])[0]) mexErrMsgTxt("Guide image must have the same height than input image.");
  /* longueur image */
  dims[1] = mxGetDimensions(prhs[0])[1];
  if (dims[1] != mxGetDimensions(prhs[1])[1]) mexErrMsgTxt("Guide image must have the same width than input image.");
  /* nombre de composantes couleurs */
  if (ndims == 3) dims[2] = mxGetDimensions(prhs[0])[2];
  else dims[2]=1;
  ndimsGuide = mxGetNumberOfDimensions(prhs[1]);
  if (ndimsGuide == 3) ncol=mxGetDimensions(prhs[1])[2];
  else ncol=1;

  /* allocation of the result image */
  plhs[0] = mxCreateNumericArray(ndims, dims, mxUINT8_CLASS, mxREAL);
  /* pointer on the image result buffer */
  result = (unsigned char *)mxGetData(plhs[0]);

  /* bilateral filter algorithm */
  switch (dims[2]) {
	case 1:	/* gray level */
 		if (!JointCrossBilateralFilter(dims[0], dims[1], ncol, input, guide, hwsize, sscale, gscale, result))
			mexErrMsgTxt("Error in BilateralFilter.");    
		break;
	case 3: /* color */
		/* red */
 		if (!JointCrossBilateralFilter(dims[0], dims[1], ncol, input, guide, hwsize, sscale, gscale, result))
			mexErrMsgTxt("Error in BilateralFilter.");    
		/* green */
 		if (!JointCrossBilateralFilter(dims[0], dims[1], ncol, input+dims[0]*dims[1], guide, hwsize, sscale, 
			gscale, result+dims[0]*dims[1])) mexErrMsgTxt("Error in BilateralFilter.");    
		/* blue */
 		if (!JointCrossBilateralFilter(dims[0], dims[1], ncol, input+2*dims[0]*dims[1], guide, hwsize, sscale, 
			gscale, result+2*dims[0]*dims[1])) mexErrMsgTxt("Error in BilateralFilter.");  
		break;
	default:
		mexErrMsgTxt("Input image should be in gray level or RGB.");
 }

}

