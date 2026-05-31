function Result = collimation_analysis(FitsFile, Args)
% COLLIMATION_ANALYSIS  Collimation analysis for an annular-aperture telescope.
%
%  Reads an out-of-focus star-field FITS image and evaluates telescope
%  collimation by analyzing the annular (donut) star images.
%
%  REQUIREMENTS
%    - MATLAB Image Processing Toolbox
%    - fitsread / fitsinfo  (built into MATLAB)
%
%  USAGE
%    % Run with all defaults:
%    Result = collimation_analysis('/path/to/image.fits');
%
%    % Override selected parameters:
%    Args.PixelScale    = 1.25;
%    Args.DetectThresh  = 8.0;
%    Result = collimation_analysis('/path/to/image.fits', Args);
%
%  INPUT
%    FitsFile  - Full path to the FITS file (char / string)
%    Args      - (optional) struct with any subset of the fields below.
%                Missing fields are filled with the listed default values.
%
%  ARGS FIELDS  (all optional – defaults shown)
%    Detector & optics
%      Args.PixelScale          = 1.25    % arcsec / pixel
%      Args.ObstructionRatio    = 0.40    % secondary / primary diameter
%    Star detection
%      Args.BgSigmaClip         = 3.0    % sigma-clipping half-width for sky bg
%      Args.DetectThresh        = 10.0   % detection threshold (sigma above bg)
%      Args.StarRadiusMin       = 8      % minimum expected star radius (px)
%      Args.StarRadiusMax       = 150     % maximum expected star radius (px)
%      Args.MaxStars            = 500    % analyse at most this many brightest stars
%    Ellipse fitting
%      Args.MinRingPixels       = 60     % minimum ring-edge pixels for valid fit
%    Collimation thresholds
%      Args.OffsetWarnPix       = 1.5    % ring-centre offset (px) → WARNING
%      Args.OffsetFailPix       = 4.0    % ring-centre offset (px) → POOR
%      Args.OffsetWarnNorm      = 0.05   % normalised offset → WARNING
%      Args.OffsetFailNorm      = 0.15   % normalised offset → POOR
%      Args.EllipWarn           = 0.05   % ellipticity (1-b/a) → WARNING
%      Args.EllipFail           = 0.12   % ellipticity → POOR
%    Annulus quality filter
%      Args.MinContrastRatio    = 2.0    % ring_median / hole_median minimum
%      Args.MinAnnularSnr       = 6.0   % ring SNR minimum
%      Args.MaxHoleFillFrac     = 0.45   % hole_median / ring_median maximum
%      Args.MinRingCompleteness = 0.65   % azimuthal sector fill fraction minimum
%      Args.MinBlobCircularity  = 0.40   % blob area / bounding-circle area minimum
%   Verbose
%      Args.Verbose = false             % Script actions in detail
%      Args.DisplayFigs = false         % if true displays the figs and
%      saves a fig file as well

%  OUTPUTS  (all saved in <fits_dir>/<fits_basename>/)
%    - result_table.txt                   : per-star metrics
%    - star_quality.txt                   : per-star annulus quality filter log
%    - fig1_annotated_field.fig/.png      : annotated star field
%    - fig2_offset_map.fig/.png           : spatial offset scatter map
%    - fig3a_offset_histogram_px.fig/.png : offset distribution (pixels)
%    - fig3b_offset_histogram_norm.fig/.png : offset distribution (normalised)
%    - fig4_misalignment_vectors.fig/.png : polar misalignment plot
%    - fig5_quality_metrics.fig/.png      : annulus quality metrics (5 subplots)
%
%  RETURN VALUE
%    Result  - struct with summary statistics (means, stds, verdict, paths)
%
%  Written by Claude Sonnet 4.6, April 2026

arguments
    FitsFile
    % Detector & optics
    Args.PixelScale          = 1.25
    Args.ObstructionRatio    = 0.40
    % Star detection
    Args.BgSigmaClip         = 3.0
    Args.DetectThresh        = 30.0
    Args.StarRadiusMin       = 8
    Args.StarRadiusMax       = 150
    Args.MaxStars            = 500
    % Ellipse fitting
    Args.MinRingPixels       = 60
    % Collimation thresholds
    Args.OffsetWarnPix       = 1.5
    Args.OffsetFailPix       = 4.0
    Args.OffsetWarnNorm      = 0.05
    Args.OffsetFailNorm      = 0.15
    Args.EllipWarn           = 0.05
    Args.EllipFail           = 0.12
    % Annulus quality filter
    Args.MinContrastRatio    = 1.5
    Args.MinAnnularSnr       =  3.0
    Args.MaxHoleFillFrac     = 0.45
    Args.SectorFillThreshold = 0.5
    Args.MinRingCompleteness = 0.65
    Args.MinBlobCircularity  = 0.30
    % Verbose
    Args.Verbose             = false
    Args.DisplayFigs         = false
end
% Unpack Args into local variables for readability throughout the function
PIXEL_SCALE           = Args.PixelScale;
OBSTRUCTION_RATIO     = Args.ObstructionRatio;
BG_SIGMA_CLIP         = Args.BgSigmaClip;
DETECT_THRESH         = Args.DetectThresh;
STAR_RADIUS_MIN       = Args.StarRadiusMin;
STAR_RADIUS_MAX       = Args.StarRadiusMax;
MAX_STARS             = Args.MaxStars;
MIN_RING_PIXELS       = Args.MinRingPixels;
OFFSET_WARN_PIX       = Args.OffsetWarnPix;
OFFSET_FAIL_PIX       = Args.OffsetFailPix;
OFFSET_WARN_NORM      = Args.OffsetWarnNorm;
OFFSET_FAIL_NORM      = Args.OffsetFailNorm;
ELLIP_WARN            = Args.EllipWarn;
ELLIP_FAIL            = Args.EllipFail;
MIN_CONTRAST_RATIO    = Args.MinContrastRatio;
MIN_ANNULAR_SNR       = Args.MinAnnularSnr;
MAX_HOLE_FILL_FRAC    = Args.MaxHoleFillFrac;
SECTOR_FILL_THRESHOLD = Args.SectorFillThreshold;
MIN_RING_COMPLETENESS = Args.MinRingCompleteness;
MIN_BLOB_CIRCULARITY  = Args.MinBlobCircularity;
VERBOSE               = Args.Verbose;
DISPLAY_FIGS          = Args.DisplayFigs;


%% =========================================================================
%  1. LOAD FITS FILE
%  =========================================================================
if isempty(FitsFile)
    error('collimation_analysis: FitsFile argument is required.');
end

[~, fits_name_, fits_ext_] = fileparts(FitsFile);
if VERBOSE, fprintf('Loading: %s\n', [fits_name_ fits_ext_]); end
% info  = fitsinfo(FitsFile);
raw   = double(fitsread(FitsFile));

if ndims(raw) == 3
    raw = raw(:,:,1);
    warning('3-D array detected - using first plane only.');
end

[nrows, ncols] = size(raw);
if VERBOSE, fprintf('Image size: %d x %d pixels\n', ncols, nrows); end

%% =========================================================================
%  1b. SET UP RESULT FOLDER
%  =========================================================================
% Derive the base name (no extension) from FitsFile and create a results
% sub-folder inside the same directory as the FITS file.
[fits_dir, fits_base, ~] = fileparts(FitsFile);
RESULT_DIR = fullfile(fits_dir, fits_base);

if ~exist(RESULT_DIR, 'dir')
    mkdir(RESULT_DIR);
    if VERBOSE, fprintf('Created result folder: %s\n', RESULT_DIR); end
else
    if VERBOSE, fprintf('Result folder already exists: %s\n', RESULT_DIR); end
end

%% =========================================================================
%  2. PREPROCESSING
%  =========================================================================
% if size(raw, 3) == 3
%     raw = 0.2126*raw(:,:,1) + 0.7152*raw(:,:,2) + 0.0722*raw(:,:,3);
% end

% Sigma-clipped sky background
img_flat = raw(:);
bg_med   = median(img_flat);
bg_std   = std(img_flat);
for k = 1:5
    mask   = abs(img_flat - bg_med) < BG_SIGMA_CLIP * bg_std;
    bg_med = median(img_flat(mask));
    bg_std = std(img_flat(mask));
end
bg_std = max(bg_std, 1.0);  % floor at 1 ADU to prevent degenerate detection
if VERBOSE, fprintf('Sky background: %.1f  Noise (sigma): %.2f\n', bg_med, bg_std); end

img = raw - bg_med;

%% =========================================================================
%  3. STAR DETECTION
%  =========================================================================
thresh_map   = img > DETECT_THRESH * bg_std;
thresh_fill  = imfill(thresh_map, 'holes');
min_area     = pi * STAR_RADIUS_MIN^2;
max_area     = pi * STAR_RADIUS_MAX^2;
thresh_clean = bwareaopen(thresh_fill, round(min_area));

CC   = bwconncomp(thresh_clean);
rp   = regionprops(CC, img, ...
       'Centroid','Area','MajorAxisLength','MinorAxisLength', ...
       'Orientation','WeightedCentroid','BoundingBox','Solidity');

valid = ([rp.Area] > min_area) & ([rp.Area] < max_area);
rp    = rp(valid);

if isempty(rp)
    error(['No stars detected. Try adjusting DETECT_THRESH, ' ...
           'STAR_RADIUS_MIN/MAX, or check the FITS file.']);
end

% Sort by flux (descending)
flux = zeros(1, numel(rp));
for k = 1:numel(rp)
    bb  = round(rp(k).BoundingBox);
    r1  = max(1, bb(2)); r2 = min(nrows, bb(2)+bb(4));
    c1  = max(1, bb(1)); c2 = min(ncols, bb(1)+bb(3));
    flux(k) = sum(sum(img(r1:r2, c1:c2)));
end
[~, sidx] = sort(flux, 'descend');
rp         = rp(sidx(1:min(MAX_STARS, numel(sidx))));
nstars     = numel(rp);
if VERBOSE, fprintf('Stars selected for analysis: %d\n', nstars); end

%% =========================================================================
%  4. RING EXTRACTION, QUALITY FILTERING, AND ELLIPSE FITTING
%  =========================================================================
n_fail_contrast    = 0;
n_fail_snr         = 0;
n_fail_hole        = 0;
n_fail_complete    = 0;
n_fail_circularity = 0;

% Open star_quality.txt for writing (written inside the loop below)
sq_file = fullfile(RESULT_DIR, 'star_quality.txt');
fid_sq  = fopen(sq_file, 'w');
if fid_sq == -1
    warning('Could not open star_quality.txt for writing: %s', sq_file);
end
% cleanup_sq = onCleanup(@() fclose_if_open(fid_sq)); 
if fid_sq ~= -1
    fprintf(fid_sq, '  %-8s %-6s %-10s %-8s %-12s %-12s %-10s  %s\n', ...
            'Star', 'Status', 'Contrast', 'SNR', 'HoleFillFrac', 'Complete', 'Circ', 'FailReasons');
    fprintf(fid_sq, '  %s\n', repmat('-', 1, 75));
end
results = repmat(struct('id',[],'qc',[],'cx_outer',[],'cy_outer',[], ...
    'cx_inner',[],'cy_inner',[],'ellipticity_outer',[],'ellipticity_inner',[], ...
    'ep_outer_stamp',[],'ep_inner_stamp',[],'outer_radius',[],'offset_pix',[], ...
    'offset_norm',[],'offset_arcsec',[],'angle_offset',[],'radial_variance',[]), ...
    nstars, 1);
for k = 1:nstars
    bb  = rp(k).BoundingBox;
    pad = round(0.3 * max(bb(3), bb(4)));
    r1  = max(1, round(bb(2)) - pad);
    r2  = min(nrows, round(bb(2)+bb(4)) + pad);
    c1  = max(1, round(bb(1)) - pad);
    c2  = min(ncols, round(bb(1)+bb(3)) + pad);

    stamp    = img(r1:r2, c1:c2);
    [sr, sc] = size(stamp);

    stamp_norm   = mat2gray(stamp);
    stamp_filled = imfill(stamp_norm, 'holes');

    % Approximate geometric centre of this blob (in stamp coords)
    cx_s = rp(k).WeightedCentroid(1) - c1 + 1;
    cy_s = rp(k).WeightedCentroid(2) - r1 + 1;

    % Outer radius estimate from bounding box
    r_outer_est = 0.5 * max(bb(3), bb(4));
    r_inner_est = r_outer_est * OBSTRUCTION_RATIO;

    [Xg, Yg] = meshgrid(1:sc, 1:sr);
    Rg        = sqrt((Xg - cx_s).^2 + (Yg - cy_s).^2);

    annulus_mask = (Rg >= r_inner_est * 0.8) & (Rg <= r_outer_est * 1.1);
    hole_mask_q  = Rg < r_inner_est * 0.7;

    ring_pixels = stamp(annulus_mask);
    hole_pixels = stamp(hole_mask_q);

    ring_median  = median(ring_pixels);
    ring_mean    = mean(ring_pixels);
    if ~isempty(hole_pixels)
        hole_median  = median(hole_pixels);
    else
        hole_median = 0;
    end

    if hole_median <= 0
        contrast_ratio = Inf;
    else
        contrast_ratio = ring_median / hole_median;
    end

    annular_snr = ring_mean / (bg_std + eps);

    if ring_median > 0
        hole_fill_frac = hole_median / ring_median;
    else
        hole_fill_frac = 1.0;
    end

    n_sectors    = 16;
    sector_angle = 2*pi / n_sectors;
    theta_map    = atan2(Yg - cy_s, Xg - cx_s);
    sector_lit   = false(1, n_sectors);
    for s = 1:n_sectors
        a0 = -pi + (s-1) * sector_angle;
        a1 = a0 + sector_angle;
        sec_mask = annulus_mask & (theta_map >= a0) & (theta_map < a1);
        if any(sec_mask(:)) && median(stamp(sec_mask)) > SECTOR_FILL_THRESHOLD * ring_median
            sector_lit(s) = true;
        end
    end
    ring_completeness = sum(sector_lit) / n_sectors;

    blob_area          = rp(k).Area;
    bounding_r         = 0.5 * sqrt(bb(3)^2 + bb(4)^2);
    bounding_circ_area = pi * bounding_r^2;
    blob_circularity   = blob_area / bounding_circ_area;

    qc = struct();
    qc.contrast_ratio    = contrast_ratio;
    qc.annular_snr       = annular_snr;
    qc.hole_fill_frac    = hole_fill_frac;
    qc.ring_completeness = ring_completeness;
    qc.blob_circularity  = blob_circularity;

    fail_contrast    = contrast_ratio    < MIN_CONTRAST_RATIO;
    fail_snr         = annular_snr       < MIN_ANNULAR_SNR;
    fail_hole        = hole_fill_frac    > MAX_HOLE_FILL_FRAC;
    fail_complete    = ring_completeness < MIN_RING_COMPLETENESS;
    fail_circularity = blob_circularity  < MIN_BLOB_CIRCULARITY;

    qc.pass = ~(fail_contrast || fail_snr || fail_hole || ...
                fail_complete || fail_circularity);

    n_fail_contrast    = n_fail_contrast    + fail_contrast;
    n_fail_snr         = n_fail_snr         + fail_snr;
    n_fail_hole        = n_fail_hole        + fail_hole;
    n_fail_complete    = n_fail_complete    + fail_complete;
    n_fail_circularity = n_fail_circularity + fail_circularity;

    % Write per-star quality line to star_quality.txt
    status_str = 'PASS';
    if ~qc.pass, status_str = 'FAIL'; end
    reasons = '';
    if fail_contrast,    reasons = [reasons ' contrast'  ]; end
    if fail_snr,         reasons = [reasons ' SNR'       ]; end
    if fail_hole,        reasons = [reasons ' holefill'  ]; end
    if fail_complete,    reasons = [reasons ' incomplete']; end
    if fail_circularity, reasons = [reasons ' circ'     ]; end
    if fid_sq ~= -1
        fprintf(fid_sq, '  %-8d %-6s %-10.2f %-8.1f %-10.2f %-12.2f %-10.2f  %s\n', ...
                k, status_str, contrast_ratio, annular_snr, ...
                hole_fill_frac, ring_completeness, blob_circularity, strtrim(reasons));
    end

    edges_out  = edge(stamp_filled, 'canny', [0.05 0.15]);
    hole_diff  = stamp_filled - stamp_norm;
    hole_mask2 = hole_diff > (0.15 * max(hole_diff(:)));
    edges_in   = bwperim(hole_mask2);

    [ry_o, rx_o] = find(edges_out);
    [ry_i, rx_i] = find(edges_in);

    res_k    = struct();
    res_k.id = k;
    res_k.qc = qc;

    if numel(rx_o) >= MIN_RING_PIXELS
        ep_o = fit_ellipse_algebraic(rx_o, ry_o);
        res_k.cx_outer          = ep_o.cx + c1 - 1;
        res_k.cy_outer          = ep_o.cy + r1 - 1;
        res_k.ellipticity_outer = 1 - ep_o.minor / ep_o.major;
        res_k.ep_outer_stamp    = ep_o;
        res_k.outer_radius    = (ep_o.major + ep_o.minor) / 2;   % mean of semi-axes = mean outer radius
    else
        ep_o = struct('minor', r_inner_est, 'major', r_outer_est, ...
                      'cx', cx_s, 'cy', cy_s, 'angle', 0);
        res_k.cx_outer          = rp(k).WeightedCentroid(1);
        res_k.cy_outer          = rp(k).WeightedCentroid(2);
        res_k.ellipticity_outer = 1 - rp(k).MinorAxisLength / ...
                                       max(rp(k).MajorAxisLength, 1);
        res_k.ep_outer_stamp    = [];
        res_k.outer_radius    = (r_outer_est + r_inner_est) / 2;  % fallback estimate
        warning('Star %d: insufficient outer ring pixels for ellipse fit.', k);
    end

    if numel(rx_i) >= MIN_RING_PIXELS
        ep_i = fit_ellipse_algebraic(rx_i, ry_i);
        res_k.cx_inner          = ep_i.cx + c1 - 1;
        res_k.cy_inner          = ep_i.cy + r1 - 1;
        res_k.ellipticity_inner = 1 - ep_i.minor / ep_i.major;
        res_k.ep_inner_stamp    = ep_i;
    else
        res_k.cx_inner          = rp(k).WeightedCentroid(1);
        res_k.cy_inner          = rp(k).WeightedCentroid(2);
        res_k.ellipticity_inner = res_k.ellipticity_outer;
        res_k.ep_inner_stamp    = [];
        warning('Star %d: insufficient inner ring pixels for ellipse fit.', k);
    end

    dx = res_k.cx_outer - res_k.cx_inner;
    dy = res_k.cy_outer - res_k.cy_inner;
    res_k.offset_pix    = sqrt(dx^2 + dy^2);
    res_k.offset_norm   = res_k.offset_pix / max(res_k.outer_radius, eps); % dimensionless: offset / outer_radius
    res_k.offset_arcsec = res_k.offset_pix * PIXEL_SCALE;
    res_k.angle_offset  = atan2d(dy, dx);

    cx_s2 = res_k.cx_outer - c1 + 1;
    cy_s2 = res_k.cy_outer - r1 + 1;
    res_k.radial_variance = radial_brightness_variance(stamp, cx_s2, cy_s2, ...
                                round(ep_o.minor * 0.6), ...
                                round(ep_o.major * 1.1));

    results(k) = res_k;
end

% Close star_quality.txt
if fid_sq ~= -1
    fclose(fid_sq);
    if VERBOSE, fprintf('Star quality saved: %s\n', sq_file); end
end

%% =========================================================================
%  5. COLLIMATION ASSESSMENT  (quality-filtered stars only)
%  =========================================================================
pass_idx = find(arrayfun(@(r) r.qc.pass, results));
fail_idx = find(arrayfun(@(r) ~r.qc.pass, results));

n_pass = numel(pass_idx);
n_fail = numel(fail_idx);

fprintf('Quality filter: %d detected, %d passed, %d rejected\n', nstars, n_pass, n_fail);
if VERBOSE
    fprintf('  Fail-contrast  : %d\n',   n_fail_contrast);
    fprintf('  Fail-SNR       : %d\n',   n_fail_snr);
    fprintf('  Fail-holefill  : %d\n',   n_fail_hole);
    fprintf('  Fail-complete  : %d\n',   n_fail_complete);
    fprintf('  Fail-circularity: %d\n',  n_fail_circularity);
end

if n_pass == 0
    error(['All stars failed the quality filter. ' ...
           'Try relaxing the MIN_CONTRAST_RATIO, MIN_ANNULAR_SNR, ' ...
           'MIN_RING_COMPLETENESS or other thresholds.']);
end

offsets_pix_all    = [results.offset_pix];
offsets_norm_all   = [results.offset_norm];
ellipticities_all  = [results.ellipticity_outer];
outer_rad_all      = [results.outer_radius];

offsets_pix_pass   = offsets_pix_all(pass_idx);
offsets_norm_pass  = offsets_norm_all(pass_idx);
ellipticities_pass = ellipticities_all(pass_idx);
outer_rad_pass    = outer_rad_all(pass_idx);

mean_outer_rad     = mean(outer_rad_pass);           % mean outer_radius over passing stars (px)
std_outer_rad      = std(outer_rad_pass);            % std of outer_radius (px)
mean_offset        = mean(offsets_pix_pass);         % mean offset in pixels
std_offset         = std(offsets_pix_pass);          % std of offset (px)
median_offset      = median(offsets_pix_pass);       % median offset in pixels
mean_offset_norm   = mean(offsets_norm_pass);        % mean normalised offset (dimensionless)
std_offset_norm    = std(offsets_norm_pass);         % std of normalised offset
median_offset_norm = median(offsets_norm_pass);      % median normalised offset (dimensionless)
mean_ellip         = mean(ellipticities_pass);       % mean ellipticity
std_ellip          = std(ellipticities_pass);        % std of ellipticity

cx_outer_pass = [results(pass_idx).cx_outer];
cy_outer_pass = [results(pass_idx).cy_outer];
cx_inner_pass = [results(pass_idx).cx_inner];
cy_inner_pass = [results(pass_idx).cy_inner];

dx_pass    = cx_outer_pass - cx_inner_pass;
dy_pass    = cy_outer_pass - cy_inner_pass;
mean_dx    = mean(dx_pass);
mean_dy    = mean(dy_pass);
std_dx     = std(dx_pass);                           % std of x-component of offset
std_dy     = std(dy_pass);                           % std of y-component of offset
mean_angle = atan2d(mean_dy, mean_dx);
std_angle  = std(atan2d(dy_pass, dx_pass));          % std of per-star misalignment angle (deg)

if mean_offset_norm < OFFSET_WARN_NORM && mean_ellip < ELLIP_WARN
    verdict = 'GOOD- collimation is within acceptable limits';
elseif mean_offset_norm < OFFSET_FAIL_NORM && mean_ellip < ELLIP_FAIL
    verdict = 'WARNING- slight collimation error detected';
else
    verdict = 'POOR- significant collimation error detected';
end

%% =========================================================================
%  6. RESULTS TABLE  (console + saved text file)
%  =========================================================================
report_lines = {};   % accumulate lines for file output

hdr = sprintf('==============================================================');
report_lines{end+1} = hdr;
report_lines{end+1} = '  COLLIMATION REPORT  (quality-filtered stars only)';
report_lines{end+1} = hdr;
report_lines{end+1} = sprintf('  File  : %s', FitsFile);
report_lines{end+1} = sprintf('  Stars : %d analysed, %d passed quality filter', nstars, n_pass);
report_lines{end+1} = sprintf('  Mean outer radius       : %.2f +/- %.2f px', mean_outer_rad, std_outer_rad);
report_lines{end+1} = sprintf('  Mean ring-centre offset : %.2f +/- %.2f px  (%.2f +/- %.2f arcsec)  [%.4f +/- %.4f norm]', ...
                               mean_offset, std_offset, ...
                               mean_offset * PIXEL_SCALE, std_offset * PIXEL_SCALE, ...
                               mean_offset_norm, std_offset_norm);
report_lines{end+1} = sprintf('  Median offset           : %.2f px  [%.4f norm]', ...
                               median_offset, median_offset_norm);
report_lines{end+1} = sprintf('  Mean outer ellipticity  : %.3f +/- %.3f  (0 = perfect circle)', mean_ellip, std_ellip);
report_lines{end+1} = sprintf('  Mean misalignment angle : %.1f +/- %.1f deg (image coords)', mean_angle, std_angle);
report_lines{end+1} = sprintf('');
report_lines{end+1} = sprintf('  VERDICT: %s', verdict);
report_lines{end+1} = hdr;
report_lines{end+1} = '';
report_lines{end+1} = sprintf('  Quality Filter Summary');
report_lines{end+1} = sprintf('  Total detected  : %d', nstars);
report_lines{end+1} = sprintf('  Passed          : %d', n_pass);
report_lines{end+1} = sprintf('  Rejected        : %d', n_fail);
report_lines{end+1} = sprintf('  Fail-contrast   : %d', n_fail_contrast);
report_lines{end+1} = sprintf('  Fail-SNR        : %d', n_fail_snr);
report_lines{end+1} = sprintf('  Fail-holefill   : %d', n_fail_hole);
report_lines{end+1} = sprintf('  Fail-complete   : %d', n_fail_complete);
report_lines{end+1} = sprintf('  Fail-circularity: %d', n_fail_circularity);
report_lines{end+1} = '';

% n_summary_lines marks the end of the summary section (printed to console).
% The per-star table rows appended after this point go to the file only.
n_summary_lines = numel(report_lines);

col_hdr = sprintf('  %-5s %-5s %-12s %-10s %-10s %-12s %-14s %-10s %-12s %-10s', ...
                  'Star', 'Qual', 'OuterRad(px)', 'Offset(px)', 'OffsetNorm', 'Offset(")', ...
                  'Ellipticity', 'RadVar', 'Angle(deg)', 'Contrast');
report_lines{end+1} = col_hdr;
report_lines{end+1} = sprintf('  %s', repmat('-', 1, 100));

for k = 1:nstars
    qual_str = 'OK';
    if ~results(k).qc.pass, qual_str = '--'; end
    row = sprintf('  %-5d %-5s %-12.2f %-10.2f %-10.4f %-12.2f %-14.3f %-10.4f %-12.1f %-10.2f', ...
                  k, qual_str, results(k).outer_radius, results(k).offset_pix, ...
                  results(k).offset_norm, results(k).offset_arcsec, ...
                  results(k).ellipticity_outer, results(k).radial_variance, ...
                  results(k).angle_offset, results(k).qc.contrast_ratio);
    report_lines{end+1} = row;
end
report_lines{end+1} = '';

% Print summary to console (header + stats + quality filter; no per-star table)
if VERBOSE
    for i = 1:n_summary_lines
        fprintf('%s\n', report_lines{i});
    end
end

% Save to text file
txt_file = fullfile(RESULT_DIR, 'result_table.txt');
fid = fopen(txt_file, 'w');
if fid == -1
    warning('Could not open result table file for writing: %s', txt_file);
else
    for i = 1:numel(report_lines)
        fprintf(fid, '%s\n', report_lines{i});
    end
    fclose(fid);
    if VERBOSE, fprintf('Result table saved: %s\n', txt_file); end
end

%% =========================================================================
%  7. DIAGNOSTIC PLOTS  -  one figure per plot
%  =========================================================================
% save_fig(fig, base) -> saves fig as base.fig and base.png  [local function below]

% Common colour definitions reused across figures
BG_FIG  = [0.12 0.12 0.15];   % figure background
BG_AX   = [0.15 0.15 0.20];   % axes background

% -------------------------------------------------------------------------
% Figure 1 – Annotated star field
% -------------------------------------------------------------------------
fig1 = figure('Name', 'Annotated Star Field', 'NumberTitle', 'off', ...
              'Color', BG_FIG, 'Position', [50 50 1100 850], ...
              'Visible', bool2visstr(DISPLAY_FIGS));
ax1  = axes(fig1);
display_img = log1p(max(img, 0));
imshow(display_img, [], 'Parent', ax1);
colormap(ax1, 'gray');
hold(ax1, 'on');
title(ax1, sprintf('Annotated star field  (pass = colour-coded, fail = red)  |  Verdict: %s', verdict), ...
      'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');

for k = 1:nstars
    is_pass = results(k).qc.pass;
    if is_pass
        pass_rank = find(pass_idx == k);
        if n_pass > 1
            clr = jet(n_pass);
            clr = clr(pass_rank, :);
        else
            clr = [0.3 0.9 0.3];
        end
        lw_outer = 1.8;  lw_inner = 1.0;
    else
        clr      = [0.85 0.2 0.2];
        lw_outer = 1.0;  lw_inner = 0.7;
    end
    if ~isempty(results(k).ep_outer_stamp)
        draw_ellipse(ax1, results(k).cx_outer, results(k).cy_outer, ...
                     results(k).ep_outer_stamp.major, ...
                     results(k).ep_outer_stamp.minor, ...
                     results(k).ep_outer_stamp.angle, clr, lw_outer);
    end
    if ~isempty(results(k).ep_inner_stamp)
        draw_ellipse(ax1, results(k).cx_inner, results(k).cy_inner, ...
                     results(k).ep_inner_stamp.major, ...
                     results(k).ep_inner_stamp.minor, ...
                     results(k).ep_inner_stamp.angle, clr, lw_inner);
    end
    if is_pass
        plot(ax1, [results(k).cx_outer, results(k).cx_inner], ...
                  [results(k).cy_outer, results(k).cy_inner], ...
                  'Color', [1 0.4 0.1], 'LineWidth', 1.2);
    end
    text(ax1, results(k).cx_outer+5, results(k).cy_outer-5, ...
         sprintf('%d', k), 'Color', clr, 'FontSize', 8);
end
hold(ax1, 'off');
save_fig(fig1, fullfile(RESULT_DIR, 'fig1_annotated_field'), VERBOSE , DISPLAY_FIGS);

% -------------------------------------------------------------------------
% Figure 2 – Ring offset spatial map (scatter over detector footprint)
% -------------------------------------------------------------------------
fig2 = figure('Name', 'Ring Offset Map', 'NumberTitle', 'off', ...
              'Color', BG_FIG, 'Position', [200 100 700 600], ...
              'Visible', bool2visstr(DISPLAY_FIGS));
ax2  = axes(fig2);
set(ax2, 'Color', BG_AX, 'XColor', 'w', 'YColor', 'w');
hold(ax2, 'on');
if n_pass > 0
    scatter(ax2, cx_outer_pass, cy_outer_pass, 25, offsets_pix_pass, 'filled');
    cb = colorbar(ax2);
    cb.Color = 'w';
    cb.Label.String = 'Ring centre offset (px)';
    cb.Label.Color  = 'w';
    colormap(ax2, 'hot');
end
if n_fail > 0
    scatter(ax2, [results(fail_idx).cx_outer], [results(fail_idx).cy_outer], ...
            25, 'rx', 'LineWidth', 1.5);
    legend(ax2, 'Pass (offset)', 'Fail (rejected)', ...
           'TextColor', 'w', 'Color', BG_AX, 'EdgeColor', 'w');
end
title(ax2, sprintf('Ring offset map  (%d/%d stars pass)', n_pass, nstars), ...
      'Color', 'w', 'FontSize', 11);
xlabel(ax2, 'X (px)', 'Color', 'w');
ylabel(ax2, 'Y (px)', 'Color', 'w');
axis(ax2, 'image');
set(ax2, 'YDir', 'reverse');
hold(ax2, 'off');
save_fig(fig2, fullfile(RESULT_DIR, 'fig2_offset_map'), VERBOSE , DISPLAY_FIGS);

% -------------------------------------------------------------------------
% Figure 3a – Ring centre offset histogram (pixels)
% -------------------------------------------------------------------------
fig3a = figure('Name', 'Offset Histogram (px)', 'NumberTitle', 'off', ...
               'Color', BG_FIG, 'Position', [350 150 700 500], ...
               'Visible', bool2visstr(DISPLAY_FIGS));
ax3a  = axes(fig3a);
set(ax3a, 'Color', BG_AX, 'XColor', 'w', 'YColor', 'w');
hold(ax3a, 'on');
if n_pass > 1
    histogram(ax3a, offsets_pix_pass, max(5, round(n_pass/2)), ...
              'FaceColor', [0.3 0.6 1], 'EdgeColor', 'none');
elseif n_pass == 1
    bar(ax3a, offsets_pix_pass, 'FaceColor', [0.3 0.6 1]);
end
xline(ax3a, OFFSET_WARN_PIX, '--', 'Color', [1 0.7 0], 'LineWidth', 1.8, ...
      'Label', sprintf('Warn (%.1f px)', OFFSET_WARN_PIX), 'LabelColor', [1 0.7 0]);
xline(ax3a, OFFSET_FAIL_PIX, '--', 'Color', [1 0.3 0.1], 'LineWidth', 1.8, ...
      'Label', sprintf('Fail (%.1f px)', OFFSET_FAIL_PIX), 'LabelColor', [1 0.3 0.1]);
xline(ax3a, mean_offset, '-', 'Color', [0.4 1 0.4], 'LineWidth', 2.0, ...
      'Label', sprintf('Mean (%.2f px)', mean_offset), 'LabelColor', [0.4 1 0.4]);
xlabel(ax3a, 'Ring centre offset (px)', 'Color', 'w');
ylabel(ax3a, 'Count', 'Color', 'w');
title(ax3a, sprintf('Offset distribution  (pixels)  |  %d / %d stars used', n_pass, nstars), ...
      'Color', 'w', 'FontSize', 11);
hold(ax3a, 'off');
save_fig(fig3a, fullfile(RESULT_DIR, 'fig3a_offset_histogram_px'), VERBOSE , DISPLAY_FIGS);

% -------------------------------------------------------------------------
% Figure 3b – Ring centre offset histogram (normalised: offset / outer_radius)
% -------------------------------------------------------------------------
fig3b = figure('Name', 'Offset Histogram (normalised)', 'NumberTitle', 'off', ...
               'Color', BG_FIG, 'Position', [420 200 700 500], ...
               'Visible', bool2visstr(DISPLAY_FIGS));
ax3b  = axes(fig3b);
set(ax3b, 'Color', BG_AX, 'XColor', 'w', 'YColor', 'w');
hold(ax3b, 'on');
norm_bin_edges = linspace(0, 1, max(6, round(n_pass/2) + 1));  % fixed [0,1] bins
if n_pass > 1
    histogram(ax3b, offsets_norm_pass, norm_bin_edges, ...
              'FaceColor', [0.2 0.8 0.6], 'EdgeColor', 'none');
elseif n_pass == 1
    bar(ax3b, offsets_norm_pass, 'FaceColor', [0.2 0.8 0.6]);
end
xlim(ax3b, [0, 1]);
xline(ax3b, OFFSET_WARN_NORM, '--', 'Color', [1 0.7 0], 'LineWidth', 1.8, ...
      'Label', sprintf('Warn (%.3f)', OFFSET_WARN_NORM), 'LabelColor', [1 0.7 0]);
xline(ax3b, OFFSET_FAIL_NORM, '--', 'Color', [1 0.3 0.1], 'LineWidth', 1.8, ...
      'Label', sprintf('Fail (%.3f)', OFFSET_FAIL_NORM), 'LabelColor', [1 0.3 0.1]);
xline(ax3b, mean_offset_norm, '-', 'Color', [0.4 1 0.4], 'LineWidth', 2.0, ...
      'Label', sprintf('Mean (%.4f)', mean_offset_norm), 'LabelColor', [0.4 1 0.4]);
xlabel(ax3b, 'Normalised offset  (offset / outer\_radius)', 'Color', 'w');
ylabel(ax3b, 'Count', 'Color', 'w');
title(ax3b, sprintf('Offset distribution  (normalised)  |  %d / %d stars used', n_pass, nstars), ...
      'Color', 'w', 'FontSize', 11);
hold(ax3b, 'off');
save_fig(fig3b, fullfile(RESULT_DIR, 'fig3b_offset_histogram_norm'), VERBOSE,  DISPLAY_FIGS);

% -------------------------------------------------------------------------
% Figure 4 – Misalignment polar plot
% -------------------------------------------------------------------------
fig4 = figure('Name', 'Misalignment Vectors', 'NumberTitle', 'off', ...
              'Color', BG_FIG, 'Position', [500 200 600 560], ...
              'Visible', bool2visstr(DISPLAY_FIGS));
ax4  = polaraxes('Parent', fig4);
set(ax4, 'Color', BG_AX, 'ThetaColor', 'w', 'RColor', 'w', ...
         'GridColor', [0.5 0.5 0.5]);
hold(ax4, 'on');
% Fixed radial range 0-1 and threshold reference circles
ax4.RLim = [0, 1];
ax4.RTick = 0:0.2:1;
theta_circle = linspace(0, 2*pi, 361);
polarplot(ax4, theta_circle, repmat(OFFSET_WARN_NORM, 1, 361), ...
          '--', 'Color', [1 0.7 0], 'LineWidth', 1.2);    % warn ring
polarplot(ax4, theta_circle, repmat(OFFSET_FAIL_NORM, 1, 361), ...
          '--', 'Color', [1 0.3 0.1], 'LineWidth', 1.2);  % fail ring
cmap_pass = jet(max(n_pass, 1));
for qi = 1:n_pass
    k = pass_idx(qi);
    polarplot(ax4, deg2rad(results(k).angle_offset), results(k).offset_norm, ...
              'o', 'Color', cmap_pass(qi,:), 'MarkerSize', 4, ...
              'MarkerFaceColor', cmap_pass(qi,:));
end
for qi = 1:numel(fail_idx)
% for k = fail_idx
    k= fail_idx(qi);
    polarplot(ax4, deg2rad(results(k).angle_offset), results(k).offset_norm, ...
              'd', 'Color', [0.5 0.5 0.5], 'MarkerSize', 4);
    % polarplot(ax4, deg2rad(results(k).angle_offset), results(k).offset_norm, ...
    %           'd', 'MarkerEdgeColor', [0.5 0.5 0.5], ...
    %           'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 4);
end
polarplot(ax4, [0, deg2rad(mean_angle)], [0, mean_offset_norm], '-w', 'LineWidth', 2.5);
polarplot(ax4, deg2rad(mean_angle), mean_offset_norm, ...
          'p', 'Color', 'w', 'MarkerSize', 4, 'MarkerFaceColor', 'y');
title(ax4, sprintf('Misalignment vectors  (normalised offset)  |  %d used, %d rejected', n_pass, n_fail), ...
      'Color', 'w', 'FontSize', 11);
hold(ax4, 'off');
save_fig(fig4, fullfile(RESULT_DIR, 'fig4_misalignment_vectors'), VERBOSE, DISPLAY_FIGS);

% -------------------------------------------------------------------------
% Figure 5 – Annulus quality metrics (5 subplots in one row)
% -------------------------------------------------------------------------
qc_fields  = {'contrast_ratio', 'annular_snr', 'hole_fill_frac', ...
              'ring_completeness', 'blob_circularity'};
qc_labels  = {'Contrast Ratio', 'Annular SNR', 'Hole Fill Fraction', ...
              'Ring Completeness', 'Blob Circularity'};
qc_thresh  = [MIN_CONTRAST_RATIO, MIN_ANNULAR_SNR, MAX_HOLE_FILL_FRAC, ...
              MIN_RING_COMPLETENESS, MIN_BLOB_CIRCULARITY];
qc_invert  = [false, false, true, false, false];
qc_thresh_label = {'Min', 'Min', 'Max', 'Min', 'Min'};

nq    = numel(qc_fields);
fig5  = figure('Name', 'Annulus Quality Metrics', 'NumberTitle', 'off', ...
               'Color', BG_FIG, 'Position', [50 50 1400 380], ...
               'Visible', bool2visstr(DISPLAY_FIGS));
for qi = 1:nq
    ax_q = subplot(1, nq, qi);
    set(ax_q, 'Color', BG_AX, 'XColor', 'w', 'YColor', 'w', 'Parent', fig5);

    vals      = arrayfun(@(r) r.qc.(qc_fields{qi}), results);
    vals_disp = min(vals, 3 * qc_thresh(qi));

    bar_colors = zeros(nstars, 3);
    for k = 1:nstars
        if results(k).qc.pass
            bar_colors(k,:) = [0.3 0.6 1.0];           % blue  = all checks pass
        else
            if qc_invert(qi)
                culprit = vals(k) > qc_thresh(qi);
            else
                culprit = vals(k) < qc_thresh(qi);
            end
            if culprit
                bar_colors(k,:) = [1.0 0.3 0.2];       % red   = this metric failed
            else
                bar_colors(k,:) = [0.6 0.6 0.6];       % grey  = failed elsewhere
            end
        end
    end

    b        = bar(ax_q, 1:nstars, vals_disp, 'FaceColor', 'flat');
    b.CData  = bar_colors;
    yline(ax_q, qc_thresh(qi), '--y', 'LineWidth', 1.8, ...
          'Label', sprintf('%s = %.2f', qc_thresh_label{qi}, qc_thresh(qi)), ...
          'LabelColor', 'y');
    xlabel(ax_q, 'Star #', 'Color', 'w');
    if qi == 1, ylabel(ax_q, 'Value', 'Color', 'w'); end
    title(ax_q, qc_labels{qi}, 'Color', 'w', 'FontSize', 10);
    xlim(ax_q, [0, nstars+1]);
end
save_fig(fig5, fullfile(RESULT_DIR, 'fig5_quality_metrics'), VERBOSE, DISPLAY_FIGS);

if VERBOSE
    fprintf('\nAnalysis complete. Adjust the secondary mirror opposite to\n');
    fprintf('the mean misalignment direction (%.1f deg in image coords).\n', ...
            mean_angle);
    fprintf('Quality filter kept %d of %d stars.\n', n_pass, nstars);
    fprintf('\nAll results written to: %s\n\n', RESULT_DIR);
end

%% =========================================================================
%  POPULATE RETURN STRUCT
%  =========================================================================
Result.FitsFile          = FitsFile;
Result.ResultDir         = RESULT_DIR;
Result.NStars            = nstars;
Result.NPass             = n_pass;
Result.NFail             = n_fail;
Result.MeanOuterRad      = mean_outer_rad;
Result.StdOuterRad       = std_outer_rad;
Result.MeanOffsetPix     = mean_offset;
Result.StdOffsetPix      = std_offset;
Result.MedianOffsetPix   = median_offset;
Result.MeanOffsetNorm    = mean_offset_norm;
Result.StdOffsetNorm     = std_offset_norm;
Result.MedianOffsetNorm  = median_offset_norm;
Result.MeanOffsetArcsec  = mean_offset * PIXEL_SCALE;
Result.StdOffsetArcsec   = std_offset  * PIXEL_SCALE;
Result.MeanEllipticity   = mean_ellip;
Result.StdEllipticity    = std_ellip;
Result.MeanDx            = mean_dx;
Result.StdDx             = std_dx;
Result.MeanDy            = mean_dy;
Result.StdDy             = std_dy;
Result.MeanAngle         = mean_angle;
Result.StdAngle          = std_angle;
Result.Verdict           = verdict;

% Save Args to file in ResultDir
args_file = fullfile(RESULT_DIR, 'args.txt');
fid_args  = fopen(args_file, 'w');
if fid_args ~= -1
    fields = fieldnames(Args);
    fprintf(fid_args, 'collimation_analysis  Args\n');
    fprintf(fid_args, '%-25s  %s\n', 'Parameter', 'Value');
    fprintf(fid_args, '%s\n', repmat('-', 1, 50));
    for fi = 1:numel(fields)
        val = Args.(fields{fi});
        if islogical(val)
            valstr = mat2str(val);
        elseif isnumeric(val)
            valstr = num2str(val);
        else
            valstr = char(val);
        end
        fprintf(fid_args, '%-25s  %s\n', fields{fi}, valstr);
    end
    fclose(fid_args);
    if VERBOSE, fprintf('Args saved: %s\n', args_file); end
else
    warning('Could not write args file: %s', args_file);
end

end  % function collimation_analysis

%% =========================================================================
%  LOCAL FUNCTIONS
%  =========================================================================

function s = bool2visstr(tf)
%BOOL2VISSTR  Return 'on' or 'off' from a logical value.
    if tf, s = 'on'; else, s = 'off'; end
end

function save_fig(f, base, verbose,displayfigs)
%SAVE_FIG  Save figure F as both a .fig and a .png file.
    if displayfigs
        savefig(f, [base '.fig']);
    end
    exportgraphics(f, [base '.png'], 'Resolution', 150)
    if nargin < 3, verbose = false; end
    if verbose, fprintf('Saved: %s.fig / .png\n', base); 
    end
end

function fclose_if_open(fid)
%FCLOSE_IF_OPEN  Close a file handle only if it is still open.
%  Used by onCleanup to guarantee fid_sq is released even if the
%  function throws before reaching the explicit fclose() call.
    if fid ~= -1
        try
            fclose(fid);
        catch
            % already closed — nothing to do
        end
    end
end

function ep = fit_ellipse_algebraic(x, y)
%FIT_ELLIPSE_ALGEBRAIC  Least-squares algebraic ellipse fit (Fitzgibbon 1999)
    x = x(:); y = y(:);
    D = [x.^2, x.*y, y.^2, x, y, ones(size(x))];
    S = D' * D;
    C = zeros(6,6); C(1,3) = 2; C(3,1) = 2; C(2,2) = -1;
    [V, lam] = eig(S \ C);
    lambda   = diag(lam);
    pos_def  = (lambda > 1e-10) & isfinite(lambda);
    if ~any(pos_def)
        ep = fallback_ellipse(x, y); return;
    end
    [~, idx] = min(lambda(pos_def));
    vp = V(:, pos_def);
    a  = vp(:, idx);
    ep = conic_to_geometric(a);
end

function ep = conic_to_geometric(a)
%CONIC_TO_GEOMETRIC  Convert conic coefficients to geometric ellipse params
    A = a(1); B = a(2)/2; C = a(3); D = a(4)/2; F = a(5)/2; G = a(6);
    denom = B^2 - A*C;
    if abs(denom) < 1e-12
        ep = fallback_ellipse([], []); return;
    end
    cx  = (C*D - B*F) / denom;
    cy  = (A*F - B*D) / denom;
    num = 2*(A*F^2 + C*D^2 + G*B^2 - 2*B*D*F - A*C*G);
    s1  = sqrt((A-C)^2 + 4*B^2);
    a2  = -sqrt(num / (denom*(s1 - A - C)));
    b2  = -sqrt(num / (denom*(-s1 - A - C)));
    if isnan(a2) || isnan(b2) || ~isreal(a2) || ~isreal(b2)
        ep = fallback_ellipse([], []); return;
    end
    ep.cx    = cx;
    ep.cy    = cy;
    ep.major = max(abs(a2), abs(b2));
    ep.minor = min(abs(a2), abs(b2));
    ep.angle = 0.5 * atan2d(2*B, A-C);
end

function ep = fallback_ellipse(x, y)
%FALLBACK_ELLIPSE  Simple bounding-box estimate when algebraic fit fails
    if isempty(x)
        ep.cx = 0; ep.cy = 0; ep.major = 10; ep.minor = 10; ep.angle = 0;
    else
        ep.cx    = mean(x); ep.cy = mean(y);
        ep.major = (max(x)-min(x))/2;
        ep.minor = (max(y)-min(y))/2;
        ep.angle = 0;
    end
end

function rv = radial_brightness_variance(stamp, cx, cy, r_in, r_out)
%RADIAL_BRIGHTNESS_VARIANCE  Azimuthal variance of the bright ring
    [rows, cols] = size(stamp);
    [X, Y] = meshgrid(1:cols, 1:rows);
    R    = sqrt((X - cx).^2 + (Y - cy).^2);
    ring = stamp(R >= r_in & R <= r_out);
    if isempty(ring)
        rv = 0; return;
    end
    rv = std(ring) / (mean(ring) + eps);
end

function draw_ellipse(ax, cx, cy, major, minor, angle_deg, color, lw)
%DRAW_ELLIPSE  Overlay an ellipse on axes AX
    t   = linspace(0, 2*pi, 200);
    ang = deg2rad(angle_deg);
    xe  = major * cos(t);
    ye  = minor * sin(t);
    xr  = cx + xe * cos(ang) - ye * sin(ang);
    yr  = cy + xe * sin(ang) + ye * cos(ang);
    plot(ax, xr, yr, 'Color', color, 'LineWidth', lw);
end