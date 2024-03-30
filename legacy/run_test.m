
% read image
imafilename = 'hudson_diatomE.pgm';
%imafilename = 'peppers.ppm';
orig = imread(imafilename);

% read guide image
guidefilename = 'hudson_diatomE-Guide.pgm';
%guidefilename = 'peppers-Guide.ppm';
guide = imread(guidefilename);
psnr_guide=psnr(orig,guide)

% gaussian noise
nstd=10.0;
noisy1 = uint8(double(orig)+randn(size(orig))*nstd);
figure(1), imshow(noisy1), title('With Gaussian noise');
psnr_noisy1=psnr(orig,noisy1)

% gaussian and salt and pepper noise
nper=0.1;
noisy2 = saltpepper(noisy1,nper);
figure(2), imshow(noisy2), title('With Gaussian and salt and pepper noise');
psnr_noisy2=psnr(orig,noisy2)

% half size of the filtering window
hwsize=2;
% spatial scale
sscale=1.0;

bf_res1=bilateralfilter(noisy1, hwsize, sscale, nstd);
figure(3), imshow(bf_res1), title('Bilateral filter on Gaussian noise'); 
psnr_bf_res1=psnr(orig,bf_res1)

bf_res2=bilateralfilter(noisy2, hwsize, sscale, nstd);
figure(4), imshow(bf_res2), title('Bilateral filter on non Gaussian noise'); 
psnr_bf_res2=psnr(orig,bf_res2)

rbf_res1=robustbilateralfilter(noisy1, hwsize, sscale, nstd, 0.5);
figure(5), imshow(rbf_res1), title('Robust bilateral filter on Gaussian noise'); 
psnr_rbf_res1=psnr(orig,rbf_res1)

rbf_res2=robustbilateralfilter(noisy2, hwsize, sscale, nstd, 0.5);
figure(6), imshow(rbf_res2), title('Robust bilateral filter on non Gaussian noise'); 
psnr_rbf_res2=psnr(orig,rbf_res2)

jcbf_res1=jointcrossbilateralfilter(noisy1, guide, hwsize, sscale, nstd);
figure(7), imshow(jcbf_res1), title('Joint/cross bilateral filter on Gaussian noise'); 
psnr_jcbf_res1=psnr(orig,jcbf_res1)

jcbf_res2=jointcrossbilateralfilter(noisy2, guide, hwsize, sscale, nstd);
figure(8), imshow(jcbf_res2), title('Joint/cross bilateral filter on non Gaussian noise'); 
psnr_jcbf_res2=psnr(orig,jcbf_res2)

gbf_res1=guidedbilateralfilter(noisy1, guide, 5, 0.0, nstd, 0.5, min(10,2*nstd), 1.0);
figure(9), imshow(gbf_res1), title('Guided bilateral filter on Gaussian noise'); 
psnr_gbf_res1=psnr(orig,gbf_res1)

gbf_res2=guidedbilateralfilter(noisy2, guide, 5, 0.0, nstd, 0.5, min(10,2*nstd), 0.0);
figure(10), imshow(gbf_res2), title('Guided bilateral filter on non Gaussian noise'); 
psnr_gbf_res2=psnr(orig,gbf_res2)



