function resu=saltpepper(ima,percent)

resu=ima;
uni=rand(size(ima));
ind=find(uni<percent);
resu(ind)=255;
ind=find(uni<percent*0.5);
resu(ind)=0;



