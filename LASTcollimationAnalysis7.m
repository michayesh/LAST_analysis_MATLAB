clear;
close all;
clc;
FITS_PATH         = '/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST.01.08.04_focus_20250421/';
FITS_FILE         = fullfile(FITS_PATH,'LAST.01.08.04_20250421.165429.218_clear_138+60_001_001_001_focus_raw_Image_1.fits');

Result = collimation_analysis(FITS_FILE,'Verbose',true,'DisplayFigs',true);