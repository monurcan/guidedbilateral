/*
 * =============================================================
 * robustbilateralfilter_mex.c 
 *
 * Matlab-mex interface to the robust bilateral filtering 
 * of an image
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2015 IFSTTAR
 * 
 * Jean-Philippe TAREL (11 february 2015)
 *
 * =============================================================
 */

/* $Revision: 0.0.0.1 $ */

#include <string.h>
#include <math.h>

#include "mex.h"

int RobustBilateralFilterStep(int dimx, int dimy, unsigned char *orig, int demisize, float sscale, 
	float iscale, float ipower, float *filtered)
{
 int i, j, k, l, value, ediff;
 float somme, poids, pixelMoy, currentIntensity, diff, rdiff, *sweight=NULL, iweight[257];    

 /* spatial weight */
 if ( (sweight = (float *)malloc( (demisize+1) * sizeof(float) ) ) == NULL) return(0);
 for(i=0; i<=demisize; i++) {
	if (sscale>0.0f) sweight[i]=exp(-0.5f*(float)(i*i)/(sscale*sscale));
	else sweight[i]=1.0f;
 }

 /* intensity weight */
 for(i=0; i<=256; i++) {
	if (ipower!=1.0f) iweight[i]=pow(1.0f+(float)(i*i)/(iscale*iscale),ipower-1.0f);
	else iweight[i]=1.0f;
 }

 for(j=0; j<dimy; j++) {
 	for(i=0; i<dimx; i++) {
		somme = 1e-6f;
		pixelMoy = 0.0f;
		currentIntensity=filtered[j*dimx+i];
		for(k=-demisize; k<=demisize; k++) {  
			if ((j+k>=0) && (j+k<dimy)) {
				for(l=-demisize; l<=demisize; l++) {
			  		if ((i+l>=0) && (i+l<dimx)) {
			    			value = orig[(j+k)*dimx+i+l];
			    			diff = fabs((float)value-currentIntensity);
						ediff=(int)floor(diff);
						rdiff=diff-(float)ediff;
						poids = ((1.0f-rdiff)*iweight[ediff]+rdiff*iweight[ediff+1])*sweight[abs(k)]*sweight[abs(l)]; 
			    			somme += poids;
			    			pixelMoy += poids*(float)value;
					}
			  	}
			}
		}
		filtered[j*dimx+i]=pixelMoy/somme;
	}
 }
 free(sweight);

 return(1);
}


int RobustBilateralFilter(int dimx, int dimy, unsigned char *orig, int demisize, float sscale, float iscale, float ipower, unsigned char *result)
{
 int i, num=8;
 float *filtered=NULL;

 /* alloc */
 if ( (filtered = malloc( dimx * dimy * sizeof(float) ) ) == NULL) return(0);

 /* init image */
 for(i=0; i<dimx*dimy; i++) filtered[i]=(float)(orig[i]);

 /* GNC */
 if (ipower<=1.0f) {
	if (!RobustBilateralFilterStep(dimx, dimy, orig, demisize, 0.0, iscale, 1.0f, filtered)) return(0);
	num--;
 }

 if (ipower<=0.5f) {
	if (!RobustBilateralFilterStep(dimx, dimy, orig, demisize, sscale, iscale, 0.5f, filtered)) return(0);
	num--;
 }
 
 if (ipower<=0.0f) {
	if (!RobustBilateralFilterStep(dimx, dimy, orig, demisize, sscale, iscale, 0.0f, filtered)) return(0);
	num--;
 }

 /* final */
 for(i=0; i<num; i++) {
	if (!RobustBilateralFilterStep(dimx, dimy, orig, demisize, sscale, iscale, ipower, filtered)) return(0);
 } 

 for(i=0; i<dimx*dimy; i++) result[i]=(unsigned char)(filtered[i]);
			
 free(filtered);

 return(1);
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  unsigned char *input=NULL, *result=NULL;
  int ndims, dims[3], hwsize=2; 
  float sscale=1.5f, iscale=10.0f, ipower=0.0f;

  /*  Check for proper number of arguments. */
  if (nrhs < 1) mexErrMsgTxt("One input argument is at least required.");

  /* Check to make sure the first input argument is an uint8 image */
  if (!mxIsClass(prhs[0],"uint8")) mexErrMsgTxt("Input image should be uint8.");

  /* Create a pointer to the input image */
  input = (unsigned char *)mxGetData(prhs[0]);

  if (nrhs > 1) hwsize=(int)mxGetScalar(prhs[1]);
  if (nrhs > 2) sscale=(float)mxGetScalar(prhs[2]);
  if (nrhs > 3) iscale=(float)mxGetScalar(prhs[3]);
  if (nrhs > 4) ipower=(float)mxGetScalar(prhs[4]);
  if (nrhs > 5) mexErrMsgTxt("Too many input arguments.");
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

  /* robust bilateral filter algorithm */
  switch (dims[2]) {
	case 1:	/* gray level */
 		if (!RobustBilateralFilter(dims[0], dims[1], input, hwsize, sscale, iscale, ipower, result))
			mexErrMsgTxt("Error in BilateralFilter.");    
		break;
	case 3: /* color */
		/* red */
 		if (!RobustBilateralFilter(dims[0], dims[1], input, hwsize, sscale, iscale, ipower, result))
			mexErrMsgTxt("Error in BilateralFilter.");    
		/* green */
 		if (!RobustBilateralFilter(dims[0], dims[1], input+dims[0]*dims[1], hwsize, sscale, iscale, ipower, result+dims[0]*dims[1]))
			mexErrMsgTxt("Error in BilateralFilter.");    
		/* blue */
 		if (!RobustBilateralFilter(dims[0], dims[1], input+2*dims[0]*dims[1], hwsize, sscale, iscale, ipower, result+2*dims[0]*dims[1]))
			mexErrMsgTxt("Error in BilateralFilter.");  
		break;
	default:
		mexErrMsgTxt("Input image should be in gray level or RGB.");
 }

}

