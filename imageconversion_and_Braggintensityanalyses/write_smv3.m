function [imageout_raw,sumframe,allpeaklabels,rezpeaklabels,peakintprofiles,radprofiles,radDvals,beamcenter]=write_smv3(datain,dose,fileinfo,scalefactor,sizeout,scalar,pedestal,SNRcutoff,write,display)
% Niko Vlahakis, Arden Clauss, and Jose Rodriguez
% NOTE this requires the sample header file (sample_smv_header.mat) to work!
% datain can be 3d matrix of Nx2D patterns.
% filename is prefix for output .img files.
% function is called by mrc2smv and tvips2smv scripts, writes out smv
% format image stacks, and measures Bragg peak intensities as a function of
% accumulated fluence

% Modification by NWV 20231101: script applies a gaussian filter only for
% peak intensity analysis, but writes out images without the filter as smv
% for processing


%% READ IN DATA FILE
%datain=double(imread(filename));
V=size(datain)*scalefactor;
halfA = floor(V(1)/2);
halfB = floor(V(2)/2);
if max(size(V)>2)
    maxframes=V(3)/scalefactor;
else
    maxframes=1;
end

Backgroundval=pedestal + 0.1;


%% USE EXISTING IMAGE TEMPLATE TO CLONE or USE TEMPLATE HEADER

% user: insert path on file system to sample_smv_header.mat in first argument here
load('#path_to_sample_smv_header.mat#','header');

%% PREPARE OUTPUT DIRECTORIES
filein  = fileinfo.path_to_mrc{1};
[filepath, filename, ~] = fileparts(filein);

%path on file system to where smv images will be written, user may change
impath = strcat(filepath,'/../images/',filename,'/');

if isfolder(impath)==0
    status = mkdir(impath);
end

%path on file system to where analyses and MATLAB variables will be stored
%for future reference, user may change
figspath=strcat(filepath,'/../analysis/',filename,'/');

if isfolder(figspath)==0
    status = mkdir(figspath);
end

%% PREPARE NEW MATRIX FOR DATA OUTPUT
if sizeout>V(1) && sizeout>V(2)
    Y=sizeout; X=sizeout;
else
    Y=max(V(1:2)); X=max(V(1:2));
end
halfY=floor(Y/2);
halfX=floor(X/2);
imageout_raw = ones(Y,X,maxframes,'single')+pedestal;


%% DEFINE EXPERIMENT PARAMETERS FOR HEADER & UPDATE TEMPLATE HEADER

%These values are for populating the smv file image metadata, and should be
%adjusted by the user to match the parameters for their experiment

pixel_size=0.0080;             % detector pix size in mm 5um for DE Apollo
wavelength=double(0.0251);     % 0.02511 A at 200KeV
osc_range=0.5;                 % currently these are filler values, can be modified if user desires accurate oscillation range in image metadata

distance=fileinfo.camera_length; % reads camera length in mm from fileinfo table

if max(Y,X)<1024
    binning=2; % compression of the image
    beam_center=[53.9 53.9]; % beam location on image, 30.7211 @ 2048
else
    binning=1; % compression of the image
    beam_center=[99.9 99.9]; % beam location on image
end
pixel_size=pixel_size*(1/scalefactor); % pixel size in mm, 0.01501 @ 4092
image_size=[Y X]; % image size

[newheader]=set_image_params(header,pixel_size,distance,osc_range,wavelength,beam_center,image_size);

for kk=1:maxframes
    %% EMBED INPUT DATA INTO OUTPUT MATRIX

    %imageout_proc is the image array modified with a Gaussian filter to
    %aid peak picking by the pick_electronPeaks function. imageout_raw is
    %the image array without this modification, which is written to smv.

    resizedimg_proc=single(imgaussfilt(imresize(datain(:,:,kk),scalefactor),1))*scalar+pedestal;
    resizedimg_raw=single(imresize(datain(:,:,kk),scalefactor))*scalar;
    if max(size(resizedimg_raw))<max(size(imageout_raw(:,:,kk)))
        imageout_proc(halfY-halfA+1:halfY+halfA,halfX-halfB+1:halfX+halfB,kk)=resizedimg_proc;
        imageout_raw(halfY-halfA+1:halfY+halfA,halfX-halfB+1:halfX+halfB,kk)=resizedimg_raw;
    else
        imageout_proc(:,:,kk) = resizedimg_proc;
        imageout_raw(:,:,kk) = resizedimg_raw;
    end
    
    if (display==1 && mod(kk,50)==0)
        [imagepeaks]=pick_electronPeaks(imageout_proc(:,:,kk),SNRcutoff,Backgroundval); %imageout
    end

    %% WRITE SMV IMG FILES
    if write==1
        if kk<10
            fid=fopen(strcat(impath,filename,'_00',num2str(kk),'.img'),'Wb');
        elseif kk<100
            fid=fopen(strcat(impath,filename,'_0',num2str(kk),'.img'),'Wb');
        else
            fid=fopen(strcat(impath,filename,'_',num2str(kk),'.img'),'Wb');
        end
        fwrite(fid,newheader,'uchar');
        %fwrite(fid,single(imageout_proc(:,:,kk)),'uint16',0,'ieee-le'); %for writing out filtered images. Comment out normally
        fwrite(fid,single(imageout_raw(:,:,kk)),'uint16',0,'ieee-le'); % imageout
        fclose(fid);
    end

    %% DISPLAY PROGRESS AS MAX DIFFRACTION PROJECTION
    imax=max(max(imageout_raw(:,:,kk)));

    if (display==1 && mod(kk,50)==0)
        figure(002),
        colormap gray;
        subplot(1,3,1), imagesc(datain(:,:,kk)), caxis([0 0.1*imax]), axis image, title(strcat('diff image: ',num2str(kk)));
        subplot(1,3,2), imagesc(imagepeaks), axis image, caxis([0 1]), title 'peaks';
        subplot(1,3,3), imagesc(imageout_raw(:,:,kk)), caxis([0 0.5*imax]), axis image, title 'diff as bin';
    end
end

%% CALCULATE MAX FRAME

%note that the variable name sumframe is outdated. This is now a MAXIMUM
%intensity projection of the provided image stack, and is used for finding
%Bragg peaks
sumframe=max(imageout_proc,[],3);

%% ANALYZE DETECTED PEAKS AND CREATE SUMMARY FIGURE
[imagepeaks,allpeaklabels,rezpeaklabels,peakintprofiles,radprofiles,radDvals,beamcenter]=radpeakanalysis2(imageout_proc,imageout_raw,dose,fileinfo,SNRcutoff,sumframe,Backgroundval);
Ang = char(197);
disp('Writing Figure 3...')
h3=figure(003); colormap jet;
subplot(2,2,1), imagesc(sumframe), caxis([0 max(max(sumframe))]), axis image, title 'Max projection of all patterns';%, set(gca,'xtick',[]);
subplot(2,2,2), imagesc(imagepeaks), axis image, caxis([0 1]), title(strcat('Total reflections: ',num2str(size(peakintprofiles,1))));%, set(gca,'xtick',[]), set(gca,'ytick',[]);
subplot(2,2,3), imagesc(allpeaklabels), caxis([0 max(max(allpeaklabels))]), axis image, title 'Labeled reflections';%, set(gca,'xtick',[]), set(gca,'ytick',[]);
subplot(2,2,4), imagesc([0 dose(end)],[],log(abs(peakintprofiles)+1)), caxis([0 max(max(log(abs(peakintprofiles)+1)))]), title 'Reflection intensity profiles over time', xlabel(sprintf('Accumulated Fluence (e-/%c^{2})',Ang)), ylabel('Peak');
saveas(h3,strcat(figspath,filename,'_f3.pdf'),'pdf');
saveas(h3,strcat(figspath,filename,'_f3.png'),'png');

%% WRITE MAX FRAME
fid=fopen(strcat(figspath,filename,'MAXFRAME.img'),'Wb');
fwrite(fid,newheader,'uchar');
fwrite(fid,single(sumframe*10/maxframes),'uint16',0,'ieee-le'); % imageout
fclose(fid);

end


%% FUNCTION TO UPDATE TEMPLATE HEADER WITH EXPT PARAMETERS
function [header]=set_image_params(header,pixel_size,distance,osc_range,wavelength,beam_center,image_size)

header_info1 = textscan(header,'%s','delimiter',';');
header_info = header_info1{1};

for kk=1:size(header_info,1)
	if isempty(header_info{kk})==0
        head_entry = textscan(header_info{kk},'%s','delimiter','=');
        temp=head_entry{1}{1};
        if strcmp(temp,'PIXEL_SIZE')==1
            head_entry{1}{2}=strcat('=',num2str(pixel_size),';');
        elseif strcmp(temp,'DISTANCE')==1
            head_entry{1}{2}=strcat('=',num2str(distance),';');
        elseif strcmp(temp,'OSC_RANGE')==1
            head_entry{1}{2}=strcat('=',num2str(osc_range),';');
        elseif strcmp(temp,'WAVELENGTH')==1
            head_entry{1}{2}=strcat('=',num2str(wavelength),';');
        elseif strcmp(temp,'BEAM_CENTER_X')==1
            head_entry{1}{2}=strcat('=',num2str(beam_center(2)),';');
        elseif strcmp(temp,'BEAM_CENTER_Y')==1
            head_entry{1}{2}=strcat('=',num2str(beam_center(1)),';');        
        elseif strcmp(temp,'CCD_IMAGE_SATURATION')==1
            head_entry{1}{2}=strcat('=',head_entry{1}{2},';');
            %sizeflag=1;
        elseif strcmp(temp,'SIZE1')==1 %&& sizeflag==1
            head_entry{1}{2}=strcat('=',num2str(image_size(1)),';');
        elseif strcmp(temp,'SIZE2')==1 %&& sizeflag==1
            head_entry{1}{2}=strcat('=',num2str(image_size(2)),';');
        elseif size(head_entry{1},1)>1
            head_entry{1}{2}=strcat('=',head_entry{1}{2},';');
        end
        header_info{kk}=char(head_entry{1});
	end
end

temp2 = char(header_info);

for kk=1:size(temp2,1)
    if kk==1
        newheader=strcat(temp2(1,:));
    end
    if mod(kk,2)==0 && kk<size(temp2,1)
        thisline = strcat(temp2(kk,:),temp2(kk+1,:));
        newheader=[newheader sprintf('%s\n',thisline)]; %strcat(newheader,temp2(kk,:),'\r\ns');
    elseif mod(kk,2)==0
        newheader=[newheader strcat(temp2(kk,:))];
    end
end

header(22:438)=newheader(21:437);

end