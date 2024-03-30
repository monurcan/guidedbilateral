function val=psnr(ima,orig)

diff=double(orig)-double(ima);
eqm=mean(mean(mean(diff.*diff)));
val=10.0*log10(255.0*255.0/eqm);



