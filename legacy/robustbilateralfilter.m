function result = robustbilateralfilter(varargin)
%ROBUSTBILATERALFILTER smooth an image while preserving strong
%    edges in presence of non-gaussian noise
%
%    ROBUSTBILATERALFILTER(I) returns the filtered image with
%    same M-by-N size as I.	
% 
%    ROBUSTBILATERALFILTER(I,hwsize) returns the filtered image
%    using a window of half size hwsize. Default value for
%    hwsize is set to 2.
%
%    ROBUSTBILATERALFILTER(I,hwsize,sscale) returns the filtered
%    image using a Gaussian spatial weight with scale sscale.
%    The weight is one when sscale is set to 0. Default value
%    for sscale is set to 1.5.
%
%    ROBUSTBILATERALFILTER(I,hwsize,sscale,iscale) returns the
%    filtered image using a SEF intensity weight with scale
%    iscale. Default value for iscale is set 10.
%
%    ROBUSTBILATERALFILTER(I,hwsize,sscale,iscale,ipower) 
%    returns the filtered image using a SEF intensity weight
%    with power ipower. Default value for ipower is set 0.
%
%    See also BILATERALFILTER, JOINTCROSSBILATERALFILTER,
%    GUIDEDBILATERALFILTER.
%
%    References:
%      [5] "Modeling Non-Gaussian Noise for Robust Image Analysis",
%      S.-S. Ieng, J.-P. Tarel, and P. Charbonnier, 
%      Proc. of International Conference on Computer Vision Theory
%      and Applications, 2007, p. 183-190.
%      
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
narginchk(1,5);

% Check validity of the input parameters 
I = varargin{1};
validateattributes(I,{'uint8'},{'real', 'nonempty', 'finite', 'nonsparse'}, mfilename,'I or RGB',1);

if (nargin>1) 
	hwsize = int32(varargin{2});
	if (numel(hwsize)~=1) || (max(hwsize(:))>min(size(I,1),size(I,2))/4) || (min(hwsize(:))<0)
		msg1 = sprintf('%s: half windows size has to be', upper(mfilename));
  		msg2 = 'a non-negative number between 0 and 1/4th of the image size.';
  		eid = sprintf('%s:outOfRangeHalfWindowSize',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	hwsize=int32(2);
end
if (nargin>2) 
	sscale=single(varargin{3});
	if (numel(sscale)~=1) || (min(sscale(:))<0)
		msg1 = sprintf('%s: spatial scale has to be', upper(mfilename));
  		msg2 = 'a non-negative number.';
  		eid = sprintf('%s:outOfRangeSpatialScale',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	sscale=single(1.5);
end
if (nargin>3) 
	iscale=single(varargin{4});
	if (numel(iscale)~=1) || (min(iscale(:))<=0)
		msg1 = sprintf('%s: intensity scale has to be', upper(mfilename));
  		msg2 = 'a positive number.';
  		eid = sprintf('%s:outOfRangeIntensityScale',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	iscale=single(10.0);
end
if (nargin>4) 
	ipower=single(varargin{5});
	if (numel(ipower)~=1) || (min(ipower(:))<=-5)
		msg1 = sprintf('%s: intensity power has to be', upper(mfilename));
  		msg2 = 'higher than -5.';
  		eid = sprintf('%s:outOfRangeIntensityPower',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	ipower=single(0.0);
end

result = robustbilateralfilter_mex(I, hwsize, sscale, iscale, ipower);

if (nargout==0) 
	imshow(result);
end


