%% COLLIMATION BATCH RUNNER
%  Calls collimation_analysis() on every FITS file in a directory,
%  collects the Result structs, and saves a CSV summary table.
%
%  USAGE
%    Set INPUT_DIR below and run the script.
%    Optionally edit the Args struct to override collimation_analysis defaults.
%
%  OUTPUT
%    <INPUT_DIR>/collimation_batch_results.csv
%
%  Written by Claude Sonnet 4.6, April 2026

clear; clc;

%% =========================================================================
%  CONFIGURATION
%  =========================================================================
INPUT_DIR = '/home/micha/Dropbox/WAO/LAST_Observations/OutOfFocusImages/LAST.01.01.02_20260430/LAST.01.01.02_20260430_oof';

% Options forwarded to collimation_analysis (edit as needed)
% Args.PixelScale          = 1.25;
% Args.ObstructionRatio    = 0.40;
% Args.DetectThresh        = 10.0;
% Args.StarRadiusMin       = 8;
% Args.StarRadiusMax       = 100;
% Args.Verbose             = false;   % suppress per-file console detail
% Args.DisplayFigs         = false;   % do not show figures on screen

%% =========================================================================
%  DISCOVER FITS FILES
%  =========================================================================
extensions = {'*.fits', '*.fit', '*.fts'};
file_list  = [];
for ei = 1:numel(extensions)
    file_list = [file_list; dir(fullfile(INPUT_DIR, extensions{ei}))]; %#ok<AGROW>
end

if isempty(file_list)
    error('No FITS files found in: %s', INPUT_DIR);
end
n_files = numel(file_list);
fprintf('Found %d FITS file(s) in %s\n\n', n_files, INPUT_DIR);

%% =========================================================================
%  LOOP OVER FILES
%  =========================================================================
% Numeric fields to collect from each Result struct (preserves column order)
num_fields = { ...
    'NStars',          'NPass',            'NFail', ...
    'MeanOuterRad',    'StdOuterRad', ...
    'MeanOffsetPix',   'StdOffsetPix',     'MedianOffsetPix', ...
    'MeanOffsetNorm',  'StdOffsetNorm',    'MedianOffsetNorm', ...
    'MeanOffsetArcsec','StdOffsetArcsec', ...
    'MeanEllipticity', 'StdEllipticity', ...
    'MeanDx',          'StdDx', ...
    'MeanDy',          'StdDy', ...
    'MeanAngle',       'StdAngle' };

n_cols      = numel(num_fields);
file_names  = cell(n_files, 1);
verdicts    = cell(n_files, 1);
num_data    = NaN(n_files, n_cols);
n_ok        = 0;
n_err       = 0;

for fi = 1:n_files
    fname     = file_list(fi).name;
    fits_path = fullfile(file_list(fi).folder, fname);

    fprintf('[%d/%d] %s ...\n ', fi, n_files, fname);

    try
        R = collimation_analysis(fits_path, 'Verbose',false,'DisplayFigs',false);

        file_names{fi} = fname;
        verdicts{fi}   = R.Verdict;

        for ci = 1:n_cols
            fld = num_fields{ci};
            if isfield(R, fld)
                num_data(fi, ci) = R.(fld);
            end
        end

        n_ok = n_ok + 1;
        fprintf('OK  [%s]\n', R.Verdict);

    catch ME
        file_names{fi} = fname;
        verdicts{fi}   = 'ERROR';
        n_err = n_err + 1;
        fprintf('FAILED: %s\n', ME.message);
    end
end

fprintf('\nBatch complete: %d succeeded, %d failed.\n', n_ok, n_err);

%% =========================================================================
%  BUILD AND SAVE CSV
%  =========================================================================
% Assemble table
T = table(file_names, verdicts, 'VariableNames', {'FileName', 'Verdict'});
for ci = 1:n_cols
    T.(num_fields{ci}) = num_data(:, ci);
end

% Save
csv_path = fullfile(INPUT_DIR, 'collimation_batch_results.csv');
writetable(T, csv_path);
fprintf('CSV saved: %s\n', csv_path);