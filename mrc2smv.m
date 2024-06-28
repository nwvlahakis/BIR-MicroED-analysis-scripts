function [imageout,sumframe,allpeaklabels,peakintprofiles]=mrc2smv(fileinfo,scalefactor,scalar,pedestal,sizeout,write,display)
%Niko Vlahakis, Arden Clauss, and Jose Rodriguez
%fileinfo = a table containing information about the dataset stored in variables sample (name of sample, e.g. "biotin"), spot_size (e.g. "10", "11"),
%   exposure_time (time in seconds described by each integrated frame of data,
%   e.g. "0.5", camera_length (calibrated camera length for the experiment in mm, e.g. "540"), and path_to_mrc (full path to the dataset in .mrc format on
%   user's file system)
%scalefactor = amount to downsample image by. e.g 0.5 --> output img is half original size
%scalar = scalar multiple applied to an image. default = 1
%pedestal = scalar additive to an image. default = 1
%sizeout = expected image size after applying scalefactor
%write = 1 for writing data, 0 for not writing data
%display = 1 for showing figures, 0 for not 

filein  = fileinfo.path_to_mrc{1};
disp('Reading File...')
mrcin = mrcread(filein);
datain = mrcin.Value;
SNRcutoff=1.25; %cutoff depends on detector sensor. ie 2.5+ for tvips xf416, or 1.25+ for apollo (DE)

[filepath, filename, ~] = fileparts(filein);

% Convert Frames to calibrated flux densities in electrons/A^2s.
% Calibration below is for UCLA Talos F200C
spotsizes = [11 10 9 8 7 6]; %Spot size settings on TEM
doserates = [0.01 0.03 0.045 0.084 0.127 0.256]; %Flux density (dose rate) in parallel beam diffraction mode at each spot size
dosetable = table(spotsizes,doserates);
%

numframes = size(datain,3);
frames = [1:numframes];
spotsize = fileinfo.spot_size;
exposure_time = fileinfo.exposure_time;
seconds = frames.*exposure_time;
fluxdensity = dosetable.doserates(dosetable.spotsizes == spotsize);
dose = seconds.*fluxdensity;

%Calls function write_smv3 to perform analyses of peak intensities and
%write out .smv format files for analysis in XDS or nXDS
[imageout, sumframe,allpeaklabels,rezpeaklabels,peakintprofiles,radprofiles,radDvals,beamcenter] = write_smv3(datain,dose,fileinfo,scalefactor,sizeout,scalar,pedestal,SNRcutoff,write,display);
end