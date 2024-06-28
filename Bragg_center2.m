function [Dmatrix,rrmatrix,beamcenter]=Bragg_center2(matrixin,fileinfo,cen,supp,showimg)
% If called by write_smv3 and radpeakanalysis2, these variables are
% specified by those functions
% matrixin = 
% cen = 
% supp = 
% showing = 

if nargin>2
    Rsupport = supp(1);
    Csupport = supp(2);
else
    Rsupport = size(matrixin,1);
    Csupport = size(matrixin,2);
end

croi=5;
rdist=160; scansteps=12; boxsize=240; stepsize=4;

%% ATTEMPT TO FIND CENTER BY CENTER-OF-MASS (SIGNAL)
if nargin>1
    tempRcenter=cen(1);
    tempCcenter=cen(2);
else
    tempRcenter=floor(Rsupport/2);
    tempCcenter=floor(Csupport/2);
    tempmatrix=matrixin(tempRcenter-croi:tempRcenter+croi,tempCcenter-croi:tempCcenter+croi);
    projy=sum(tempmatrix,2);
    projx=sum(tempmatrix,1);
    csumy=sum(projy.*(1:2*croi-1));
    csumx=sum(projx.*(1:2*croi-1));
    sall=sum(sum(tempmatrix));

    tempRcenter= round((csumy/sall)+tempRcenter);
    tempCcenter= round((csumx/sall)+tempCcenter);
end

%% OPTIMIZING THE BEAM CENTER ONE MORE TIME

roi1startR=tempRcenter-rdist-ceil(boxsize/2);
roi1endR=tempRcenter-rdist+ceil(boxsize/2);
roi1startC=tempCcenter+rdist-ceil(boxsize/2);
roi1endC=tempCcenter+rdist+ceil(boxsize/2);

roi2startR=tempRcenter+rdist-ceil(boxsize/2);
roi2endR=tempRcenter+rdist+ceil(boxsize/2);
roi2startC=tempCcenter-rdist-ceil(boxsize/2);
roi2endC=tempCcenter-rdist+ceil(boxsize/2);


tempmatrix=zeros(2*scansteps+1,2*scansteps+1,'single');
for kk=-scansteps:scansteps
    ks=kk*stepsize;
    for jj=-scansteps:scansteps
        js=jj*stepsize;
        scanmat=matrixin;
        roi1=matrixin(roi1startR+ks:roi1endR+ks,roi1startC+js:roi1endC+js);
        scanmat(roi1startR+ks:roi1endR+ks,roi1startC+js:roi1endC+js)=scanmat(roi1startR+ks:roi1endR+ks,roi1startC+js:roi1endC+js)+20;
        roi2=matrixin(roi2startR+ks:roi2endR+ks,roi2startC+js:roi2endC+js);
        scanmat(roi2startR+ks:roi2endR+ks,roi2startC+js:roi2endC+js)=scanmat(roi2startR+ks:roi2endR+ks,roi2startC+js:roi2endC+js)+40;
        roi2f=flip(flip(roi2,1),2);
        tempmatrix(kk+scansteps+1,jj+scansteps+1)=sum(sum(abs(roi1-roi2f)));
        if showimg==1 
            figure(004),
            subplot(2,2,1), imagesc(roi1), axis image;
            subplot(2,2,2), imagesc(roi2f), axis image;
            subplot(2,2,3), imagesc(tempmatrix), axis image;
            subplot(2,2,4), imagesc(scanmat), axis image, title(strcat('shift k',num2str(ks),', shift j',num2str(js)));
            drawnow; %pause(0.1);
        end
    end
end

%% identifying the shift required to move the beam center to the image center
[a, b]=min(tempmatrix,[],1);
[~, Cshift]=min(a);
Rshift=b(Cshift);

Rcenter=tempRcenter+(Rshift-scansteps-1).*stepsize;
Ccenter=tempCcenter+(Cshift-scansteps-1).*stepsize;

cenMat=matrixin; cenMat(Rcenter-10:Rcenter+10,Ccenter-10:Ccenter+10)=max(max(matrixin));

[rrmatrix, Dmatrix]=setKspacegrid(size(cenMat,1),size(cenMat,2),Rcenter,Ccenter,fileinfo);

beamcenter = horzcat(Ccenter,Rcenter);

figure(005),
subplot(1,3,1), imagesc(cenMat), title('Found Center'), axis image;
subplot(1,3,2), imagesc(rrmatrix), title('Radius Matrix'), axis image;
subplot(1,3,3), imagesc(Dmatrix), title('Resolution Matrix'), axis image;

end


function [rrmatrix, Dmatrix]=setKspacegrid(Y,X,ycen,xcen,fileinfo)
%% SETS UP PARAMETERS
res_cutoff= 1.0; %in Å


%distance, pixel size, wavelength (730mm = 1310mm; 1000mm = 1840mm)

%Specify correct chip size of detector in pixels and electron wavelength here. Values shown are true for UCLA
%Talos F200C with DE Apollo detector

fullchipsize=8192;               %size of apollo detector chip    
pix= 4*(fullchipsize/Y);         
lambda=0.0251;                   %wavelength of electrons in Angstroms

ddist= fileinfo.camera_length*10^3; %detector distance in units of microns
    
%% CALCULATES RESOLUTION IN RECIPROCAL SPACE
a=1:1:Y;
b=1:1:X;
[aa,bb]=meshgrid(a,b);
rrmatrix=round(sqrt((aa-ycen).^2+(bb-xcen).^2)); clear aa bb;
Dmatrix = get_kdist(rrmatrix,ddist,pix,lambda);
Dmatrix(Dmatrix<res_cutoff)=-1;


end

%% THIS IS FOR MEASURING DISTANCES IN RECIPROCAL SPACE
function [Dmatrix]=get_kdist(rrmatrix,ddist,pix,lambda)

thetamatrix=0.5*atan((pix.*rrmatrix)./ddist);
Dmatrix=lambda./(2*sin(thetamatrix));


end