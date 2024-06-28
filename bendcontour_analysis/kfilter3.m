function[kfilter]=kfilter3(matsize)

X=matsize(1); Y=matsize(2); Z=matsize(3);

xoffset=0; yoffset=0; zoffset=floor(Z/3);
midfx=X/25; midfy=Y/25; midfz=Z/10;

halfx=floor(X/2)+1;
halfy=floor(Y/2)+1;
halfz=floor(Z/2)+1;
fwhmy=Y/midfy; fwhmx=X/midfx; fwhmz=Z/midfz;
a=1:1:X; b=1:1:Y; c=1:1:Z;
[aa,bb,cc]=meshgrid(a,b,c);
kfilter=exp( - ( ((aa-halfx-xoffset).^2)./(2*fwhmx^2) + ((bb-halfy-yoffset).^2)./(2*fwhmy^2) + ((cc-halfz-zoffset).^2)./(2*fwhmz^2)) );
kfilter=kfilter+exp( - ( ((aa-halfx+xoffset).^2)./(2*fwhmx^2) + ((bb-halfy+yoffset).^2)./(2*fwhmy^2) + ((cc-halfz+zoffset).^2)./(2*fwhmz^2)) );
kfilter=kfilter/max(max(max(kfilter)));

figure(002), 
subplot(1,3,1),imagesc(squeeze(sum(kfilter,1))), axis image;
subplot(1,3,2),imagesc(squeeze(sum(kfilter,2))), axis image;
subplot(1,3,3),imagesc(squeeze(sum(kfilter,3))), axis image;