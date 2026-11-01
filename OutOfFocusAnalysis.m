% OutOfFocusAnalysis
% Use Eran's script 'findMeasureOutFocusAnnulus_m' 
% Micha May 2026
%% Read the fits file with the raw iamge 
clear;
close all;

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
%% Preprocess image
% [Status, TableRaw, AllSI, MS, Coadd, OnlyMP, JD] = pipeline.last.pipes.pipelineI(RawImageList, CI, Args)
Args.prePrepArgs                   = {};
AI = pipeline.generic.prePrep(ai, Args.prePrepArgs{:},'GlobalBackLevel',false,'MinNim',1);  
%% Do the calibrations
Args.basicCalibArgs                = {};
% AI = pipeline.generic.basicCalib(AI, CI, Args.basicCalibArgs{:}, 'UpdateJD',false,'Flat',flatImage,'Dark',darkImage);
AI = pipeline.generic.basicCalib(AI, Args.basicCalibArgs{:}, CI,'UpdateJD',false); 
%% Divide the image into 24 non-overlapping subImages (crops)
% Remove blind pixel areas
croppedImage = AI.Image(1:9576,1:6388);
% ImageSize = size(croppedImage) ;
ImageSize = [6388 9576];
SubSize = [1597 1596];
[SubImageCCDSEC, ~, ~, ~, ~] = imUtil.cut.gridSubImage(ImageSize, SubSize);
nSubIm = size(SubImageCCDSEC,1);
%% First try findMeasureOutFocusAnnulus
R=findMeasureOutFocusAnnulus_m(AI.Image,'SubImagesCCDSEC',SubImageCCDSEC ,'MaxRad',100,'Angle',15,...
        'RingsTempBank',[11 20; 15 33; 21 50; 31 60 ;50 100],'Threshold',6,'ShowFigs',true);

%%
for isub = 1:1:nSubIm
    SubImage = imUtil.cut.trim2d(croppedImage,SubImageCCDSEC(isub,:));
    % SubImgStr = sprintf("Subimage %d",isub);
    %         fig0 = figure('Name',SubImgStr, 'NumberTitle', 'off', ...
    %                'Position', [50 50 500,1050]);
    %         climits = prctile(SubImage(:), [5, 95]);
    %         imshow(SubImage,climits);
    
    Result = collimationAnalysis(SubImage,'DetectThresh',3);
end

%% Old code
% Image = fitsread(f,'primary');
% R=imUtil.psf.findMeasureOutFocusAnnulus(Image,'MaxRad',100,'Angle',15, 'RingsTempBank',[11 20; 15 33; 21 50; 31 60 ;50 100]);
% SubImagesCCDSEC   =               [2701 3700 4301 5300;...
%                                    101 1100  101 1100;...
%                                    101 1100 8501 9500;...
%                                   5301 6300  101 1100;...
%                                   5301 6300 8501 9500];
% SubImagesCCDSEC   =               [2951 3450 4551 5050;...
%                                    251 750  251 750;...
%                                    251 750 8751 9250;...
%                                   5551 6050  101 1100;...
%                                   5551 6050 8751 9250];
% R=findMeasureOutFocusAnnulus_m(f,'MaxRad',50,'Angle',5,...
%                         'HeightFactor',0.3,...
%                         'SubImagesCCDSEC', SubImagesCCDSEC,...
%                         'ShowFigs',true);
% % plot(R(1).RadialVec, R(1).MeanSector);
% fig1  = figure('Name', 'Annulus Radii', 'NumberTitle', 'off', ...
%                'Position', [50 50 2000,400]);
% % ax1 = axes(Parent=fig1)
% for i=1:numel(R)
%     subplot(1,numel(R),i);
%     polarplot(deg2rad(R(i).AngleVec(:)),R(i).InnerRadius(:));
%     hold on;
%     polarplot(deg2rad(R(i).AngleVec(:)),R(i).OuterRadius(:));
%     hold off;
% end
% Central section
% plot(R(1).AngleVec(:), [R(1).InnerRadius(:), R(1).OuterRadius(:)])
% ylim([0,1.1*max(R(1).OuterRadius)]);
% hold off;
% plot(R(1).AngleVec(:), [R(1).InnerRadius(:), R(1).OuterRadius(:)])