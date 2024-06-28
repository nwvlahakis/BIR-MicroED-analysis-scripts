function [imagepeaks,allpeaklabels,rezpeaklabels,peakintprofiles,radprofiles,radDvals,beamcenter]=radpeakanalysis2(matrixin_proc,matrixin_raw,dose,fileinfo,SNRcutoff,sumframe,Backgroundval)
% Niko Vlahakis, Arden Clauss, and Jose Rodriguez

disp('Analyzing File...')
filein  = fileinfo.path_to_mrc{1};

dshells=floor(size(sumframe,2)/2)-1; %this intentionally downsamples num of shells by 2

%% FIND PEAKS IN SUM IMAGE
[imagepeaks]=pick_electronPeaks(sumframe,SNRcutoff,Backgroundval);

[allpeaklabels,npeaks]=bwlabel(imagepeaks);
peakintprofiles=zeros(npeaks,size(matrixin_raw,3),'single');
radprofiles=zeros(dshells,size(matrixin_raw,3),'single');
radDvals=zeros(dshells,1,'single');

for kk=1:size(matrixin_raw,3)
    currpat=matrixin_raw(:,:,kk);
    for jj=1:npeaks 
        peakintprofiles(jj,kk)=sum(currpat(allpeaklabels==jj));
    end
end

%% FIND IMAGE CENTER AND CALCULATE D-SPACINGS
showbraggimg=1;
dimxy=max(size(sumframe)); halfdimxy=floor(max(size(sumframe)+1)/2);
[Dmatrix,rrmatrix,beamcenter]=Bragg_center2(sumframe,fileinfo,[halfdimxy halfdimxy],[dimxy dimxy],showbraggimg); %dimxy

peakdspacings=zeros(npeaks,1,'single');
rezpeaklabels=allpeaklabels;
for jj=1:npeaks
    peakdspacings(jj)=mean(Dmatrix(allpeaklabels==jj));
    rezpeaklabels(allpeaklabels==jj)=Dmatrix(allpeaklabels==jj);
end


%% PERFORM STATS ON MEASURED INTENSITIES
peakstds=std(peakintprofiles,[],2);
peaksmean=mean(peakintprofiles,2);

netpeakint=sum(peakintprofiles,1);
tframe=(1:size(matrixin_raw,3));

pf=polyfit(dose,netpeakint,1);
ft=polyval(pf,dose);

%% CALCULATE RADIAL PROFILES

dshellsize=(max(max(rrmatrix))-min(min(rrmatrix)))/dshells;
for kk=1:size(matrixin_raw,3)
    currpat=matrixin_raw(:,:,kk);
    currpat(allpeaklabels<1)=0;

    for jj=1:dshells
        currdmin=min(min(rrmatrix))+(jj-1)*dshellsize;
        currdmax=min(min(rrmatrix))+jj*dshellsize;
        radprofiles(jj,kk)=sum(currpat(rrmatrix>=currdmin & rrmatrix<currdmax));
        if kk==1
            radDvals(jj)=mean(Dmatrix(rrmatrix>=currdmin & rrmatrix<currdmax));
        end
    end
end