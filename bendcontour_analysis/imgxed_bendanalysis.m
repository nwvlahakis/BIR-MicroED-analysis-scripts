function [filteredstack] = imgxed_bendanalysis(pathin,pathout,flux,numcycles,framepercycle,imgsize,diffsize,pixelsize)
% Niko Vlahakis
% Function for processing paired imaging and diffraction movies in .mrc
%format to analyze bend contour motion in static microcrystals
%FUNCTION INPUTS
%pathin: path to directory where mrc files are stored
%pathout: path to directory to write processed mrc files to
%flux: incident flux density in electrons per square Angstrom per second
%numcycles:
%framepercycle:
%imgsize: dimensions in pixels of the imaging-mode micrographs (1024, 2048,
%etc.)
%diffsize: dimensions in pixels of diffraction mode images
%pixelsize: calibrated pixel size in nm in imaging mode

pathsplit = strsplit(pathin,'/');
parentdir = char(pathsplit(end));

imgstack = zeros(imgsize,imgsize,numcycles*framepercycle);


%% Accumulate imaging frames from all cycles into the same stack

for aa = 1:numcycles
    img_mrclist = dir(strcat(pathin,'/',parentdir,'_cycle',num2str(aa),'img/*.mrc'));
    img_mrcpath = strcat(pathin,'/',parentdir,'_cycle',num2str(aa),'img/');

    imgmrc = mrcread(strcat(img_mrcpath,img_mrclist(1).name));

    imgstack(:,:,((framepercycle*(aa-1))+1):framepercycle*aa) = imgmrc.Value;
end

%% Correct for bulk sample motion between batches/cycles with normalized cross correlation

alignedstack = zeros(imgsize,imgsize,numcycles*framepercycle);

alignedstack(:,:,1:framepercycle) = imgstack(:,:,1:framepercycle);
for bb = 1:numcycles-1
    referencesum = sum(alignedstack(:,:,((framepercycle*(bb-1))+1):framepercycle*bb),3);
    movingsum = sum(imgstack(:,:,((framepercycle*bb)+1):framepercycle*(bb+1)),3);
    xcorr = normxcorr2(movingsum,referencesum);
    [max_c,imax] = max(abs(xcorr(:)));
    [ypeak,xpeak] = ind2sub(size(xcorr),imax(1));
    offset = [(xpeak-size(movingsum,2)) ypeak-size(movingsum,1)];

    for cc = 1:framepercycle
        alignedstack(:,:,((framepercycle*bb)+cc)) = imtranslate(imgstack(:,:,((framepercycle*bb)+cc)),offset,'FillValues',mean(mean(mean(imgstack))));
    end
end

unalignedstack_bin = imresize3(imgstack,0.2);
alignedstack_bin = imresize3(alignedstack,0.2);


mrcwrite(alignedstack_bin,'name',strcat(pathout,'/',parentdir,'_aligned_binned.mrc'));


%% Bin stack and apply bandpass filter along time dimension


%Pad stack with three blank frames on either end -> frames filled with the
%mean pixel value of the stack
datain_nopad = alignedstack_bin;
datain = zeros(size(alignedstack_bin,1),size(alignedstack_bin,2),size(alignedstack_bin,3)+6)+mean(mean(mean(alignedstack_bin)));
datain(:,:,4:end-3) = alignedstack_bin;

kfilter=kfilter3(size(datain));

filteredstack = ifftn(ifftshift(kfilter.*fftshift(fftn(datain))));
filteredstack=filteredstack-0.5*min(min(min(filteredstack)));

figure(001)
subplot(1,3,1),imagesc(squeeze(sum(filteredstack,1))), axis image;
subplot(1,3,2),imagesc(squeeze(sum(filteredstack,2))), axis image;
subplot(1,3,3),imagesc(squeeze(sum(filteredstack,3))), axis image;

%% Write mrc and gifs for filtered stack

mrcwrite(unalignedstack_bin,'name',strcat(pathout,'/',parentdir,'_original.mrc'));
mrcwrite(filteredstack,'name',strcat(pathout,'/',parentdir,'_f3ali.mrc'));
mrcwrite(datain,'name',strcat(pathout,'/',parentdir,'_nofiltali.mrc'));


stackinmax=0.5*max(max(max(datain)));
stackoutmax=0.5*max(max(max(filteredstack)));

playspeed=0.1;
writename_og=strcat(pathout,'/',parentdir,'_ogali.gif');
writename=strcat(pathout,'/',parentdir,'_f3ali.gif');
for jj=1:size(filteredstack,3)
    if jj==1
        imwrite(uint8((100/stackinmax)*datain(:,:,jj)),writename_og,'gif','Loopcount',inf,'DelayTime',playspeed);
        imwrite(uint8((150/stackoutmax)*filteredstack(:,:,jj)),writename,'gif','Loopcount',inf,'DelayTime',playspeed);
    end
    imwrite(uint8((100/stackinmax)*datain(:,:,jj)),writename_og,'gif','WriteMode','append','Loopcount',inf,'DelayTime',playspeed);
    imwrite(uint8((150/stackoutmax)*filteredstack(:,:,jj)),writename,'gif','WriteMode','append','Loopcount',inf,'DelayTime',playspeed);
end

%% Get diffraction frames, bin, and produce filtered image montage and diffraction montage

diffstack = zeros(diffsize,diffsize,numcycles);

for ee = 1:numcycles
    diff_mrclist = dir(strcat(pathin,'/',parentdir,'_cycle',num2str(ee),'diff/*.mrc'));
    diff_mrcpath = strcat(pathin,'/',parentdir,'_cycle',num2str(ee),'diff/');

    diffmrc = mrcread(strcat(diff_mrcpath,diff_mrclist(1).name));
    diffmrc_sum = sum(diffmrc.Value,3);

    diffstack(:,:,ee) = diffmrc_sum;
end

diffstackbin = imresize(diffstack,0.2);

diffmontage2(pathout,parentdir,filteredstack(:,:,7:end-6),datain_nopad(:,:,floorDiv(size(datain_nopad,3),2)),diffstackbin,5);



%% Manually define a profile and plot bend contour propagation
%This dialogue opens an interactive figure with a projection of the image
%stack. Click two points on it to manually define a line to plot deviations
%in intensity across as a function of time/fluence


%Exclude first and last frames
filteredstack_cropped = filteredstack(:,:,7:end-6);
[stackscans,xi,yi] = linescan3(filteredstack_cropped);


%express x-axis of 2D seismogram in units of fluence
xticklabels = 0:25*flux:25*flux*(0.1*size(filteredstack_cropped,3));
xticks = linspace(1,size(stackscans,2),numel(xticklabels));


%use known pixel size to get y-axis labels for 2D seismogram in Angstroms
ytickmax = size(stackscans,1)*pixelsize*.005-mod(size(stackscans,1)*pixelsize*.005,0.5);
yticklabels = 0:0.5:ytickmax;
yticks = linspace(1,ytickmax/(pixelsize*.005),numel(yticklabels));

%get Angstrom symbol for axis labels
Ang = char(197);


f6 = figure(006);
f6.Position =  [560   397   723   316];
subplot(1,2,1), imagesc(sum(filteredstack_cropped,3)), colormap gray, axis square;
hold on;
line(xi,yi,'Color','red','LineStyle','--','LineWidth',2);
hold off;
subplot(1,2,2), imagesc(stackscans), xlabel(sprintf('Accumulated fluence (e-/%c^{2})',Ang)), ylabel(sprintf('Distance (%c)',Ang));
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'Ytick',yticks,'YTickLabel',yticklabels);

exportgraphics(f6,strcat(pathout,'/',parentdir,'_lineproj.pdf'),'ContentType','vector');
saveas(f6,strcat(pathout,'/',parentdir,'_lineproj.png'));


%% calculate absolute value of difference from the mean

ampstackscan = abs(stackscans-mean(stackscans));

amp1d = squeeze(max(ampstackscan));
amp1d_sd = squeeze(std(ampstackscan));


%test, plotting abs value change in pixel intensity for one pixel position
amp1spot= abs(stackscans(floorDiv(size(stackscans,1),2),:)-mean(stackscans));

fluence = 2.5*flux*linspace(1,size(stackscans,2),size(stackscans,2));

f7 = figure(007);
f7.Position = [436         348        1104         329];
subplot(1,3,1), imagesc(ampstackscan), colormap gray, xlabel(sprintf('Accumulated fluence (e-/%c^{2})',Ang)), ylabel(sprintf('Distance (%c)',Ang));
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'Ytick',yticks,'YTickLabel',yticklabels);
subplot(1,3,2), plot(fluence,log(amp1d)), xlabel(sprintf('Accumulated fluence (e-/%c^{2})',Ang)),ylabel log(Amplitude), xlim([0 max(fluence)]), ylim([-4 0]);
subplot(1,3,3), plot(fluence,amp1d), xlabel(sprintf('Accumulated fluence (e-/%c^{2})',Ang)),ylabel Amplitude,xlim([0 max(fluence)]), ylim([0 0.5]);
exportgraphics(f7,strcat(pathout,'/',parentdir,'_magplots_amplitudeplot.pdf'),'ContentType','vector');
saveas(f7,strcat(pathout,'/',parentdir,'_magplots_amplitudeplot.png'));

f8 = figure(008);
f8.Position = [436         348        1104         329];
subplot(1,3,1), imagesc(ampstackscan), colormap gray, xlabel(sprintf('Accumulated fluence (e-/%c^{2})',Ang)), ylabel(sprintf('Distance (%c)',Ang));
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'Ytick',yticks,'YTickLabel',yticklabels);
subplot(1,3,2), plot(fluence,amp1d_sd), xlabel(sprintf('Accumulated fluence (e-/%c^{2})',Ang)),ylabel('SD of amplitude over all line positions'), xlim([0 max(fluence)]), ylim([0 0.05]);
subplot(1,3,3), plot(fluence,amp1spot), xlabel(sprintf('Accumulated fluence (e-/%c^{2})',Ang)),ylabel('Amplitude at midpoint of line'), xlim([0 max(fluence)]), ylim([0 0.5]);
exportgraphics(f8,strcat(pathout,'/',parentdir,'_magplots_sdplot.pdf'),'ContentType','vector');
saveas(f8, strcat(pathout,'/',parentdir,'_magplots_sdplot.png'));


end