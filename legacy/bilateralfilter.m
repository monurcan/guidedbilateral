function result = bilateralfilter(varargin)
%BILATERALFILTER smooth an image while preserving strong edges
%
%    BILATERALFILTER(I) returns the filtered image with same 
%    M-by-N size as I.	
% 
%    BILATERALFILTER(I,hwsize) returns the filtered image using
%    a window of half size hwsize. Default value for hwsize is
%    set to 2.
%
%    BILATERALFILTER(I,hwsize,sscale) returns the filtered image
%    using a Gaussian spatial weight with scale sscale. The 
%    weight is one when sscale is set to 0. Default value for
%    sscale is set to 1.5.
%
%    BILATERALFILTER(I,hwsize,sscale,iscale) returns the
%    filtered image using a Gaussian intensity weight with scale
%    iscale. Default value for iscale is set 10.
%
%    See also ROBUSTBILATERALFILTER, JOINTCROSSBILATERALFILTER, 
%    GUIDEDBILATERALFILTER.
%
%    References:
%      [2] "Bilateral Filtering for Gray and Color Images",
%      C. Tomasi and R. Manduchi,
%      in ICCV, 1998, p. 839-846.
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
narginchk(1,4);

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

result = bilateralfilter_mex(I, hwsize, sscale, iscale);

if (nargout==0) 
	imshow(result);
end


