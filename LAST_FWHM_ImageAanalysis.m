

clear;
close all;
FITS_PATH  = '';
ImageName  = '';
dark_file  = ['/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST.01.04.04_20260430/calib/' ...
              'LAST.01.04.03_20241223.230027.415_clear_____dark_proc_Image_1.fits'];
flat_file  = ['/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST.01.04.04_20260430/calib/' ...
              'LAST.01.04.03_20241204.150643.935_clear_____twflat_proc_Image_1.fits'];

%% =========================================================================
%  1. READ IMAGE and calib images 
%% =========================================================================
f  = fullfile(FITS_PATH, ImageName);
ai = AstroImage(f);
fprintf('Opening %s\n', ImageName);

flatImage = AstroImage(flat_file);
darkImage = AstroImage(dark_file);
CI        = CalibImages;
CI.Dark   = darkImage;
CI.Flat   = flatImage;
CI.Bias   = darkImage;

%% =========================================================================
%  2. PREPROCESS (LAST pipeline)
%% =========================================================================
Args.prePrepArgs = {};
% Do not check PSF 
AI = pipeline.generic.prePrep(ai, Args.prePrepArgs{:}, ...
                               'GlobalBackLevel', false, 'MinNim', 1, 'GlobalBadPSF', false);

%% =========================================================================
%  3. CALIBRATE (LAST pipeline)
%% =========================================================================
Args.basicCalibArgs = {};
AI = pipeline.generic.basicCalib(AI, Args.basicCalibArgs{:}, CI, 'UpdateJD', false);

%% =========================================================================
%  4. DIVIDE INTO SUB-IMAGES (LAST pipeline)
%% =========================================================================
croppedImage = AI.Image(1:9576, 1:6388);
ImageSize    = [6388 9576];
SubSize      = [1597 1596];
[SubImageCCDSEC, ~, ~, ~, ~] = imUtil.cut.gridSubImage(ImageSize, SubSize);
nSubIm       = size(SubImageCCDSEC, 1);

