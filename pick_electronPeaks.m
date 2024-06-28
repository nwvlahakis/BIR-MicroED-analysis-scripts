function [TruePeaks]=pick_electronPeaks(mypeaksimg,SNRcutoff,Backgroundval)
% If called by write_smv3, these variables are specified by that function
% mypeaksimg = 
% SNRcutoff = 
% Backgroundval = 

Lpeaksimg=imgaussfilt(mypeaksimg,6);
Speaksimg=imgaussfilt(mypeaksimg,0.5);
Ratimg=Speaksimg./Lpeaksimg;
bordersize=20;

maskout=Speaksimg; =
maskout(maskout<Backgroundval)=0;
maskout(maskout>0)=1;
se=strel('diamond',25);
maskout=imerode(maskout,se,4);
maskout(:,1:bordersize)=0;maskout(:,size(maskout,2)-bordersize:size(maskout,2))=0;
maskout(1:bordersize,:)=0;maskout(size(maskout,1)-bordersize:size(maskout,1),:)=0;

TruePeaks=mypeaksimg;
se=strel('diamond',3);
TruePeaks(Ratimg<SNRcutoff)=0;
TruePeaks(TruePeaks>0)=1;
TruePeaks=imdilate(TruePeaks,se,2).*maskout;

figure(001),
subplot(1,3,1), imagesc(maskout), axis image, title 'Background mask'; %Lboxint
subplot(1,3,2), imagesc(Ratimg), axis image, title 'Peak targets'; %Sboxint
subplot(1,3,3), imagesc(TruePeaks.*maskout), axis image, title 'Identified peaks'; %BGboxint
