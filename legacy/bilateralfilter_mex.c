/*
 * =============================================================
 * bilateralfilter_mex.c 
 *
 * Matlab-mex interface to the bilateral filtering of an image
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2015 IFSTTAR
 * 
 * Jean-Philippe TAREL (10 february 2015)
 *
 * =============================================================
 */

/* $Revision: 0.0.0.1 $ */

#include <string.h>
#include <math.h>

#include "mex.h"

int BilateralFilter(int dimx, int dimy, unsigned char *orig, int demisize, float sscale, 
	float iscale, unsigned char *result)
{
 int i, j, k, l, value, diff, currentIntensity;
 float somme, poids, pixelMoy, *sweight=NULL, iweight[256];    

 /* spatial weight */
 if ( (sweight = (float *)malloc( (demisize+1) * sizeof(float) ) ) == NULL) return(0);
 for(i=0; i<=demisize; i++) {
	if (sscale>0.0f) sweight[i]=exp(-0.5f*((float)i*(float)i)/(sscale*sscale));
	else sweight[i]=1.0f;
 }

 /* intensity weight */
 for(i=0; i<=255; i++) iweight[i]=exp(-0.5f*(float)(i*i)/(iscale*iscale));

 for(j=0; j<dimy; j++) {
 	for(i=0; i<dimx; i++) {
		somme = 0.0f;
		pixelMoy = 0.0f;
		currentIntensity=orig[j*dimx+i];
		for(k=-demisize; k<=demisize; k++) {  
			if ((j+k>=0) && (j+k<dimy)) {
				for(l=-demisize; l<=demisize; l++) {
			  		if ((i+l>=0) && (i+l<dimx)) {
			    			value = orig[(j+k)*dimx+i+l];
			    			diff = abs(value-currentIntensity);
						poids = iweight[diff]*sweight[abs(k)]*sweight[abs(l)];
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
  unsigned char *input=NULL, *result=NULL;
  int ndims, dims[3], hwsize=2;
  float sscale=1.5f, iscale=10.0f;

  /*  Check for proper number of arguments. */
  if (nrhs < 1) mexErrMsgTxt("One input argument is at least required.");

  /* Check to make sure the first input argument is an uint8 image */
  if (!mxIsClass(prhs[0],"uint8")) mexErrMsgTxt("Input image should be uint8.");

  /* Create a pointer to the input image */
  input = (unsigned char *)mxGetData(prhs[0]);

  if (nrhs > 1) hwsize=(int)mxGetScalar(prhs[1]);
  if (nrhs > 2) sscale=(float)mxGetScalar(prhs[2]);
  if (nrhs > 3) iscale=(float)mxGetScalar(prhs[3]);
  if (nrhs > 4) mexErrMsgTxt("Too many input arguments.");
  if (nlhs < 1) mexErrMsgTxt("One output argument required.");
  if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");

  /* Get the dimensions of the image input */
  ndims = mxGetNumberOfDimensions(prhs[0]);
  /* hauteur image */
  dims[0] = mxGetDimensions(prhs[0])[0];
  /* longueur image */
  dims[1] = mxGetDimensions(prhs[0])[1];
  /* nombre de composantes couleurs */
  if (ndims == 3) dims[2] = mxGetDimensions(prhs[0])[2];
  else dims[2]=1;
 
  /* allocation of the result image */
  plhs[0] = mxCreateNumericArray(ndims, dims, mxUINT8_CLASS, mxREAL);
  /* pointer on the image result buffer */
  result = (unsigned char *)mxGetData(plhs[0]);

  /* bilateral filter algorithm */
  switch (dims[2]) {
	case 1:	/* gray level */
 		if (!BilateralFilter(dims[0], dims[1], input, hwsize, sscale, iscale, result))
			mexErrMsgTxt("Error in BilateralFilter.");    
		break;
	case 3: /* color */
		/* red */
 		if (!BilateralFilter(dims[0], dims[1], input, hwsize, sscale, iscale, result))
			mexErrMsgTxt("Error in BilateralFilter.");    
		/* green */
 		if (!BilateralFilter(dims[0], dims[1], input+dims[0]*dims[1], hwsize, sscale, iscale, result+dims[0]*dims[1]))
			mexErrMsgTxt("Error in BilateralFilter.");    
		/* blue */
 		if (!BilateralFilter(dims[0], dims[1], input+2*dims[0]*dims[1], hwsize, sscale, iscale, result+2*dims[0]*dims[1]))
			mexErrMsgTxt("Error in BilateralFilter.");  
		break;
	default:
		mexErrMsgTxt("Input image should be in gray level or RGB.");
 }

}

