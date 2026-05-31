clear;
close all;
clc;
FITS_PATH         = '/home/micha/Dropbox/WAO/LAST_Observations/OutOfFocusImages/LAST.01.04.04_20260430/failed';
[fname, fpath] = uigetfile(fullfile(FITS_PATH,'*.fits'), ...
                               'Select out-of-focus FITS image');
if isequal(fname, 0)
        error('No file selected. Aborting.');
end
    FITS_FILE = fullfile(fpath, fname);
% Result = collimation_analysis(FITS_FILE,'Verbose',true,'DisplayFigs',true);
Result = collimation_analysis(FITS_FILE,'Verbose',true,'DisplayFigs',true);