% run measure_wholeImageFWHM
% /home/micha/Dropbox/MATLABdpbx/LAST_analysis_matlab
imagespath = '/home/micha/Dropbox/WAO/LAST_analysis/psf_map/LAST01.01.03_oof';
filename = 'LAST.01.01.03_20260909.185450.241_clear_oof_36486_001_001_001_sci_raw_Image_1.fits';
[R, Aux] = measure_wholeImageFWHM_from_fits(fullfile(imagespath,filename), Method="all", ...
    Grid= [8 12],SatLevel=30000, DetectNsigma= 30);