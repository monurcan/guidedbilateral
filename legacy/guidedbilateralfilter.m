function result = guidedbilateralfilter(varargin)
%GUIDEDBILATERALFILTER smooth an image while preserving strong
%    edges in presence of non-gaussian noise with the help of a
%    guide image
%
%    GUIDEDBILATERALFILTER(I,G) returns the filtered image from
%    I with the help of guide image G. I, G and the results are 
%    with the same M-by-N size.	
% 
%    GUIDEDBILATERALFILTER(I,G,hwsize) returns the filtered
%    image using a window of half size hwsize. Default value for
%    hwsize is set to 2.
%
%    GUIDEDBILATERALFILTER(I,G,hwsize,sscale) returns the
%    filtered image using a Gaussian spatial weight with scale
%    sscale. The weight is one when sscale is set to 0. Default
%    value for sscale is set to 1.5.
%
%    GUIDEDBILATERALFILTER(I,G,hwsize,sscale,iscale) returns the
%    filtered image using a SEF intensity weight with scale
%    iscale. Default value for iscale is set 10.
%
%    GUIDEDBILATERALFILTER(I,G,hwsize,sscale,iscale,ipower) 
%    returns the filtered image using a SEF intensity weight
%    with power ipower. Default value for ipower is set 0.
%
%    GUIDEDBILATERALFILTER(I,G,hwsize,sscale,iscale,ipower,
%    gscale)
%    returns the filtered image using a SEF guide weight
%    with scale gscale. Default value for gscale is set 10.
%
%    GUIDEDBILATERALFILTER(I,G,hwsize,sscale,iscale,ipower,
%    gscale,gpower) 
%    returns the filtered image using a SEF guide weight
%    with power gpower. Default value for gpower is set 1.0.
%
%    See also BILATERALFILTER, JOINTCROSSBILATERALFILTER,
%    ROBUSTBILATERALFILTER.
%
%   References:
%      [1] "The Guided Bilateral Filter: When the Joint/Cross
%      Bilateral Filter Becomes Robust", 
%      L. Caraffa, J.-P. Tarel and P. Charbonnier,
%      in IEEE Transaction on Image Processing, 24:(4), 
%      p. 1199-1208, April 2015. 
%      http://perso.lcpc.fr/tarel.jean-philippe/publis/ip15.html
%
%   Copyright 2015 IFSTTAR.
%   $Revision: 0.0.0.2 $  $Date: 2015/07/10 14:38:00 $

% Check number of parameters 
narginchk(1,8);

% Check validity of the input parameters 
I = varargin{1};
validateattributes(I,{'uint8'},{'real', 'nonempty', 'finite', 'nonsparse'}, mfilename,'I or RGB',1);

G = varargin{2};
validateattributes(G,{'uint8'},{'real', 'nonempty', 'finite', 'nonsparse'}, mfilename,'I or RGB',2);

if (nargin>2) 
	hwsize = int32(varargin{3});
	if (numel(hwsize)~=1) || (max(hwsize(:))>min(size(I,1),size(I,2))/4) || (min(hwsize(:))<0)
		msg1 = sprintf('%s: half windows size has to be', upper(mfilename));
  		msg2 = 'a non-negative number between 0 and 1/4th of the image size.';
  		eid = sprintf('%s:outOfRangeHalfWindowSize',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	hwsize=int32(2);
end
if (nargin>3) 
	sscale=single(varargin{4});
	if (numel(sscale)~=1) || (min(sscale(:))<0)
		msg1 = sprintf('%s: spatial scale has to be', upper(mfilename));
  		msg2 = 'a non-negative number.';
  		eid = sprintf('%s:outOfRangeSpatialScale',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	sscale=single(1.5);
end
if (nargin>4) 
	iscale=single(varargin{5});
	if (numel(iscale)~=1) || (min(iscale(:))<=0)
		msg1 = sprintf('%s: intensity scale has to be', upper(mfilename));
  		msg2 = 'a positive number.';
  		eid = sprintf('%s:outOfRangeIntensityScale',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	iscale=single(10.0);
end
if (nargin>5) 
	ipower=single(varargin{6});
	if (numel(ipower)~=1) || (min(ipower(:))<=-5)
		msg1 = sprintf('%s: intensity power has to be', upper(mfilename));
  		msg2 = 'higher than -5.';
  		eid = sprintf('%s:outOfRangeIntensityPower',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	ipower=single(0.0);
end
if (nargin>5) 
	gscale=single(varargin{7});
	if (numel(gscale)~=1) || (min(gscale(:))<=0)
		msg1 = sprintf('%s: guide scale has to be', upper(mfilename));
  		msg2 = 'a positive number.';
  		eid = sprintf('%s:outOfRangeGuideScale',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	gscale=single(10.0);
end
if (nargin>7) 
	gpower=single(varargin{8});
	if (numel(gpower)~=1) || (min(gpower(:))<=-5)
		msg1 = sprintf('%s: guide power has to be', upper(mfilename));
  		msg2 = 'higher than -5.';
  		eid = sprintf('%s:outOfRangeIntensityPower',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	gpower=single(1.0);
end

result = guidedbilateralfilter_mex(I, G, hwsize, sscale, iscale, ipower, gscale, gpower);

if (nargout==0) 
	imshow(result);
end


