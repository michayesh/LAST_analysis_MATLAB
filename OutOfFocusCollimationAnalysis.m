% OutOfFocusCollimationAnalysis
% 1. Read the fits image file
% 2. Preprocess it using pipeline.generic.prePrep
% 3. Debias flatten the image using pipeline.generic.basicCalib and generic
% dark and flat images
% 4. Divide the image into 24 non-overlapping subImages (crops)
% 5. Loop over the crops and use imUtil.sources.findRings to detect the star out of focus annuli images
% 6. Choose the brightest annulus and anayze its properties 
% Micha JUne 2026
%% 1. Read the image
clear;
close all;
% Script Prameters
THRESHOLD = 6.;
RINGS_TEMP_BANK     = [11 20; 15 33; 21 50; 31 60 ; 50 100];
MAXRAD = 100;
FITS_PATH = '/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST.01.04.04_20260430/oof/';
ImageName = 'LAST.01.04.04_20260430.210246.702_clear_oof1000_001_001_001_sci_raw_Image_1.fits';
dark_file = ['/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST.01.04.04_20260430/calib/'...
   'LAST.01.04.03_20241223.230027.415_clear_____dark_proc_Image_1.fits' ];
flat_file = ['/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST.01.04.04_20260430/calib/'...
   'LAST.01.04.03_20241204.150643.935_clear_____twflat_proc_Image_1.fits' ]; 
f = fullfile(FITS_PATH,ImageName);
ai=AstroImage(fullfile(f));
fprintf('Opening %s \n',ImageName);
flatImage = AstroImage(flat_file);
darkImage = AstroImage(dark_file);
CI                      = CalibImages;
CI.Dark   = darkImage;
CI.Flat   = flatImage;
CI.Bias = darkImage;
%% 2. Preprocess image
Args.prePrepArgs                   = {};
AI = pipeline.generic.prePrep(ai, Args.prePrepArgs{:},'GlobalBackLevel',false,'MinNim',1);  
%% 3. Do the calibrations: Remove dark and bias and flatten the image
Args.basicCalibArgs                = {};
% AI = pipeline.generic.basicCalib(AI, CI, Args.basicCalibArgs{:}, 'UpdateJD',false,'Flat',flatImage,'Dark',darkImage);
AI = pipeline.generic.basicCalib(AI, Args.basicCalibArgs{:}, CI,'UpdateJD',false); 
%% 4. Divide the image into 24 non-overlapping subImages (crops)
% Remove blind pixels in the image
croppedImage = AI.Image(1:9576,1:6388);
ImageSize = size(croppedImage);
SubSize = [1597 1596];
[SubImageCCDSEC, ~, ~, ~, ~] = imUtil.cut.gridSubImage(ImageSize, SubSize);
nSubIm = size(SubImageCCDSEC,1);

%% Loop over the crops and use imUtil.sources.findRings to detect the star out of focus annuli images
for isub = 1:1:nSubIm
    SubImage = imUtil.cut.trim2d(croppedImage,SubImageCCDSEC(isub,:));
    [SizeX, SizeY]   = size(SubImage);
    [Pos,FI, Back] = ...            
            imUtil.sources.findRings(SubImage, 'CalcBack',true, 'Threshold', THRESHOLD,...
                'FunBackVar', @imUtil.background.modeVar_SampleHist,'RingsTempBank',RINGS_TEMP_BANK);
    % remove ring positions near edge:
    Flag = Pos(:,1)>MAXRAD & Pos(:,2)>MAXRAD & Pos(:,1)<(SizeX-MAXRAD) & Pos(:,2)<(SizeY-MAXRAD);
    Pos  = Pos(Flag,:);
    % Get a percentile threshold
    RingThreshold = prctile(Pos(:,3),50);
    RingPosI = find(Pos(:,3)>=RingThreshold);
    fprintf('Found %d rings in subImage %d \n',numel(RingPosI),isub);
    % Find the brightest ring in the SubImage
    [~,MaxI] = max(Pos(:,3));
    Pos_Brightest= [Pos(MaxI,1) Pos(MaxI,2)];
    Xmin = max(1, floor(Pos_Brightest(1) - MAXRAD));
    Xmax = min(SizeX, ceil(Pos_Brightest(1) + MAXRAD));
    Ymin = max(1, floor(Pos_Brightest(2) - MAXRAD));
    Ymax = min(SizeY, ceil(Pos_Brightest(2) + MAXRAD));
    % Generate the annulus cutout
    annulusCutout = double(SubImage(Ymin:Ymax, Xmin:Xmax));
end

