function result = jointcrossbilateralfilter(varargin)
%JOINTCROSSBILATERALFILTER smooth an image while preserving 
%    strong edges with the help of a guide image
%
%    JOINTCROSSBILATERALFILTER(I,G) returns the filtered image 
%    from I with the help of guide image G. I, G and the 
%    results are with the same M-by-N size.	
% 
%    JOINTCROSSBILATERALFILTER(I,G,hwsize) returns the filtered 
%    image using a window of half size hwsize. Default value for 
%    hwsize is set to 2.
%
%    JOINTCROSSBILATERALFILTER(I,G,hwsize,sscale) returns the 
%    filtered image using a Gaussian spatial weight with scale 
%    sscale. The weight is one when sscale is set to 0. Default 
%    value for sscale is set to 1.5.
%
%    JOINTCROSSBILATERALFILTER(I,G,hwsize,sscale,gscale) returns 
%    the filtered image using a Gaussian intensity guide weight 
%    with scale gscale. Default value for gscale is set 10.
%
%    See also BILATERALFILTER, ROBUSTBILATERALFILTER, 
%    GUIDEDBILATERALFILTER.
%
%    References:
%      [3] "Flash photography enhancement via intrinsic relighting",
%      E. Eisemann and F. Durand,
%      ACM Trans. on Graphics, 2004, vol. 23, num. 3, p. 673-678.
%
%      [4] "Digital photography with flash and no-flash image pairs",
%      G. Petschnigg, R. Szeliski, M. Agrawala, M. Cohen, H. Hoppe 
%      and K. Toyama, 
%      ACM Trans. on Graphics, 2004, vol. 23, num. 3, p. 664-672.
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
	gscale=single(varargin{5});
	if (numel(gscale)~=1) || (min(gscale(:))<=0)
		msg1 = sprintf('%s: intensity scale has to be', upper(mfilename));
  		msg2 = 'a positive number.';
  		eid = sprintf('%s:outOfRangeGuideScale',mfilename);
  		error(eid,'%s %s',msg1,msg2);
	end
else
	gscale=single(10.0);
end

result = jointcrossbilateralfilter_mex(I, G, hwsize, sscale, gscale);

if (nargout==0) 
	imshow(result);
end


