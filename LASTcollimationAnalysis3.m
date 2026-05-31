%% COLLIMATION ANALYSIS FOR ANNULAR APERTURE TELESCOPE
%  Reads an out-of-focus star-field FITS image and evaluates telescope
%  collimation by analyzing the annular (donut) star images.
%
%  REQUIREMENTS
%    - MATLAB Image Processing Toolbox
%    - MATLAB Astronomy & Astrophysics Toolbox  OR  a custom FITS reader
%      (the script uses fitsread / fitsinfo, which are built into MATLAB)
%
%  USAGE
%    Run the script. A file-selection dialog opens unless you set FITS_FILE
%    directly. Adjust the PARAMETERS section for your optical setup.
%
%  OUTPUTS  (all saved in <image_folder>/<fits_basename>/)
%    - result_table.txt         : per-star collimation metrics (text)
%    - collimation_analysis.fig : annotated field + offset plots (MATLAB fig)
%    - collimation_analysis.png : same figure as PNG
%    - quality_metrics.fig      : per-star annulus quality bar charts
%    - quality_metrics.png      : same figure as PNG
%    - misalignment_vectors.fig : polar misalignment plot
%    - misalignment_vectors.png : same figure as PNG
%
%  Written by Claude Sonnet 4.6, April 2026

clear all;
clc;
close all;

%% =========================================================================
%  PARAMETERS  -  edit these for your telescope and camera
%  =========================================================================
% Leave '' to open a file-selection dialog
FITS_PATH         = '/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST_20260412';
% FITS_FILE         = fullfile(FITS_PATH,'LAST.01.01.03_20260412.173436.666_clear_GoodFocusTest_001_001_001_sci_raw_Image_1.fits');
FITS_FILE         = fullfile(FITS_PATH,'LAST.01.08.02_20260412.180339.202_clear_GoodFocusTest_001_001_001_sci_raw_Image_1.fits');
% Detector & optics
PIXEL_SCALE       = 1.25;      % arcsec / pixel (for your plate scale)
OBSTRUCTION_RATIO = 0.40;      % Secondary mirror diameter / Primary diameter

% Star detection
BG_SIGMA_CLIP     = 3.0;       % Sigma-clipping iterations for background
DETECT_THRESH     = 10.0;      % Detection threshold in sigma above background
STAR_RADIUS_MIN   = 8;         % Minimum expected star radius (pixels)
STAR_RADIUS_MAX   = 50;        % Maximum expected star radius (pixels)
MAX_STARS         = 500;       % Analyse at most this many brightest stars

% Ellipse fitting
MIN_RING_PIXELS   = 60;        % Minimum pixels on each ring for a valid fit

% Collimation thresholds
OFFSET_WARN_PIX   = 1.5;       % Ring center offset (px) triggering WARNING
OFFSET_FAIL_PIX   = 4.0;       % Ring center offset (px) triggering POOR
ELLIP_WARN        = 0.05;      % Ellipticity excess (1 - b/a) for WARNING
ELLIP_FAIL        = 0.12;      % Ellipticity excess for POOR

%% =========================================================================
%  ANNULUS QUALITY FILTER PARAMETERS
%  =========================================================================
MIN_CONTRAST_RATIO    = 2.0;
MIN_ANNULAR_SNR       = 15.0;
MAX_HOLE_FILL_FRAC    = 0.45;
MIN_RING_COMPLETENESS = 0.65;
MIN_BLOB_CIRCULARITY  = 0.40;
VERBOSE_QUALITY       = true;

%% =========================================================================
%  1. LOAD FITS FILE
%  =========================================================================
if isempty(FITS_FILE)
    [fname, fpath] = uigetfile({'*.fits;*.fit;*.fts', 'FITS Files'}, ...
                               'Select out-of-focus star-field FITS image');
    if isequal(fname, 0)
        error('No file selected. Aborting.');
    end
    FITS_FILE = fullfile(fpath, fname);
end

fprintf('Loading: %s\n', FITS_FILE);
info  = fitsinfo(FITS_FILE);
raw   = double(fitsread(FITS_FILE));

if ndims(raw) == 3
    raw = raw(:,:,1);
    warning('3-D array detected - using first plane only.');
end

[nrows, ncols] = size(raw);
fprintf('Image size: %d x %d pixels\n', ncols, nrows);

%% =========================================================================
%  1b. SET UP RESULT FOLDER
%  =========================================================================
% Derive the base name (no extension) from FITS_FILE and create a results
% sub-folder inside the same directory as the FITS file.
[fits_dir, fits_base, ~] = fileparts(FITS_FILE);
RESULT_DIR = fullfile(fits_dir, fits_base);

if ~exist(RESULT_DIR, 'dir')
    mkdir(RESULT_DIR);
    fprintf('Created result folder: %s\n', RESULT_DIR);
else
    fprintf('Result folder already exists: %s\n', RESULT_DIR);
end

%% =========================================================================
%  2. PREPROCESSING
%  =========================================================================
if size(raw, 3) == 3
    raw = 0.2126*raw(:,:,1) + 0.7152*raw(:,:,2) + 0.0722*raw(:,:,3);
end

% Sigma-clipped sky background
img_flat = raw(:);
bg_med   = median(img_flat);
bg_std   = std(img_flat);
for k = 1:5
    mask   = abs(img_flat - bg_med) < BG_SIGMA_CLIP * bg_std;
    bg_med = median(img_flat(mask));
    bg_std = std(img_flat(mask));
end
fprintf('Sky background: %.1f  Noise (sigma): %.2f\n', bg_med, bg_std);

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
fprintf('Stars selected for analysis: %d\n', nstars);

%% =========================================================================
%  4. RING EXTRACTION, QUALITY FILTERING, AND ELLIPSE FITTING
%  =========================================================================
n_fail_contrast    = 0;
n_fail_snr         = 0;
n_fail_hole        = 0;
n_fail_complete    = 0;
n_fail_circularity = 0;

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
    hole_median  = median(hole_pixels);

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
        if any(sec_mask(:)) && median(stamp(sec_mask)) > 0.5 * ring_median
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

    if VERBOSE_QUALITY
        status_str = 'PASS';
        if ~qc.pass, status_str = 'FAIL'; end
        reasons = '';
        if fail_contrast,    reasons = [reasons ' contrast'   ]; end
        if fail_snr,         reasons = [reasons ' SNR'        ]; end
        if fail_hole,        reasons = [reasons ' holefill'   ]; end
        if fail_complete,    reasons = [reasons ' incomplete' ]; end
        if fail_circularity, reasons = [reasons ' shape'      ]; end
        fprintf('  Star %3d: [%s]  contrast=%.2f  SNR=%.1f holefill=%.2f  complete=%.2f  circ=%.2f %s\n', ...
                k, status_str, contrast_ratio, annular_snr, ...
                hole_fill_frac, ring_completeness, blob_circularity, reasons);
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
    else
        ep_o = struct('minor', r_inner_est, 'major', r_outer_est, ...
                      'cx', cx_s, 'cy', cy_s, 'angle', 0);
        res_k.cx_outer          = rp(k).WeightedCentroid(1);
        res_k.cy_outer          = rp(k).WeightedCentroid(2);
        res_k.ellipticity_outer = 1 - rp(k).MinorAxisLength / ...
                                       max(rp(k).MajorAxisLength, 1);
        res_k.ep_outer_stamp    = [];
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
    res_k.offset_arcsec = res_k.offset_pix * PIXEL_SCALE;
    res_k.angle_offset  = atan2d(dy, dx);

    cx_s2 = res_k.cx_outer - c1 + 1;
    cy_s2 = res_k.cy_outer - r1 + 1;
    res_k.radial_variance = radial_brightness_variance(stamp, cx_s2, cy_s2, ...
                                round(ep_o.minor * 0.6), ...
                                round(ep_o.major * 1.1));

    results(k) = res_k;
end

%% =========================================================================
%  5. COLLIMATION ASSESSMENT  (quality-filtered stars only)
%  =========================================================================
pass_idx = find(arrayfun(@(r) r.qc.pass, results));
fail_idx = find(arrayfun(@(r) ~r.qc.pass, results));

n_pass = numel(pass_idx);
n_fail = numel(fail_idx);

fprintf('\n--- Quality Filter Summary ---\n');
fprintf('  Total detected : %d\n',   nstars);
fprintf('  Passed         : %d\n',   n_pass);
fprintf('  Rejected       : %d\n',   n_fail);
fprintf('  Fail-contrast  : %d\n',   n_fail_contrast);
fprintf('  Fail-SNR       : %d\n',   n_fail_snr);
fprintf('  Fail-holefill  : %d\n',   n_fail_hole);
fprintf('  Fail-complete  : %d\n',   n_fail_complete);
fprintf('  Fail-circularity: %d\n',  n_fail_circularity);

if n_pass == 0
    error(['All stars failed the quality filter. ' ...
           'Try relaxing the MIN_CONTRAST_RATIO, MIN_ANNULAR_SNR, ' ...
           'MIN_RING_COMPLETENESS or other thresholds.']);
end

offsets_pix_all    = [results.offset_pix];
ellipticities_all  = [results.ellipticity_outer];

offsets_pix_pass   = offsets_pix_all(pass_idx);
ellipticities_pass = ellipticities_all(pass_idx);

mean_offset   = mean(offsets_pix_pass);
median_offset = median(offsets_pix_pass);
mean_ellip    = mean(ellipticities_pass);

cx_outer_pass = [results(pass_idx).cx_outer];
cy_outer_pass = [results(pass_idx).cy_outer];
cx_inner_pass = [results(pass_idx).cx_inner];
cy_inner_pass = [results(pass_idx).cy_inner];

mean_dx    = mean(cx_outer_pass - cx_inner_pass);
mean_dy    = mean(cy_outer_pass - cy_inner_pass);
mean_angle = atan2d(mean_dy, mean_dx);

if mean_offset < OFFSET_WARN_PIX && mean_ellip < ELLIP_WARN
    verdict = 'GOOD   - collimation is within acceptable limits';
    v_color = [0.0 0.6 0.2];
elseif mean_offset < OFFSET_FAIL_PIX && mean_ellip < ELLIP_FAIL
    verdict = 'WARNING - slight collimation error detected';
    v_color = [0.8 0.5 0.0];
else
    verdict = 'POOR   - significant collimation error; adjustment recommended';
    v_color = [0.8 0.1 0.1];
end

%% =========================================================================
%  6. RESULTS TABLE  (console + saved text file)
%  =========================================================================
report_lines = {};   % accumulate lines for file output

hdr = sprintf('==============================================================');
report_lines{end+1} = hdr;
report_lines{end+1} = '  COLLIMATION REPORT  (quality-filtered stars only)';
report_lines{end+1} = hdr;
report_lines{end+1} = sprintf('  File  : %s', FITS_FILE);
report_lines{end+1} = sprintf('  Stars : %d analysed, %d passed quality filter', nstars, n_pass);
report_lines{end+1} = sprintf('  Mean ring-centre offset : %.2f px  (%.2f arcsec)', ...
                               mean_offset, mean_offset * PIXEL_SCALE);
report_lines{end+1} = sprintf('  Median offset           : %.2f px', median_offset);
report_lines{end+1} = sprintf('  Mean outer ellipticity  : %.3f  (0 = perfect circle)', mean_ellip);
report_lines{end+1} = sprintf('  Mean misalignment angle : %.1f deg (image coords)', mean_angle);
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

col_hdr = sprintf('  %-5s %-5s %-10s %-12s %-14s %-10s %-12s %-10s', ...
                  'Star', 'Qual', 'Offset(px)', 'Offset(")', ...
                  'Ellipticity', 'RadVar', 'Angle(deg)', 'Contrast');
report_lines{end+1} = col_hdr;
report_lines{end+1} = sprintf('  %s', repmat('-', 1, 80));

for k = 1:nstars
    qual_str = 'OK';
    if ~results(k).qc.pass, qual_str = '--'; end
    row = sprintf('  %-5d %-5s %-10.2f %-12.2f %-14.3f %-10.4f %-12.1f %-10.2f', ...
                  k, qual_str, results(k).offset_pix, results(k).offset_arcsec, ...
                  results(k).ellipticity_outer, results(k).radial_variance, ...
                  results(k).angle_offset, results(k).qc.contrast_ratio);
    report_lines{end+1} = row;
end
report_lines{end+1} = '';

% Print to console
for i = 1:numel(report_lines)
    fprintf('%s\n', report_lines{i});
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
    fprintf('Result table saved: %s\n', txt_file);
end

%% =========================================================================
%  7. DIAGNOSTIC PLOTS
%  =========================================================================

% --- Figure 1: Collimation Analysis (annotated field + scatter + histogram)
fig1 = figure('Name', 'Collimation Analysis', 'NumberTitle', 'off', ...
       'Color', [0.12 0.12 0.15], 'Position', [50 50 1400 900]);

ax1 = subplot(2,3,[1,2,4,5]);
display_img = log1p(max(img,0));
imshow(display_img, [], 'Parent', ax1);
colormap(ax1, 'gray');
hold(ax1, 'on');
title(ax1, 'Annotated star field  (green = pass, red = fail quality)', ...
      'Color', 'w', 'FontSize', 11);

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
        lw_outer = 1.8;
        lw_inner = 1.0;
        mk_color = clr;
    else
        clr      = [0.85 0.2 0.2];
        lw_outer = 1.0;
        lw_inner = 0.7;
        mk_color = clr;
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
         sprintf('%d', k), 'Color', mk_color, 'FontSize', 8);
end
hold(ax1, 'off');

ax2 = subplot(2,3,3);
if n_pass > 0
    sc_h = scatter(cx_outer_pass, cy_outer_pass, 60, offsets_pix_pass, 'filled');
    colorbar;
    colormap(ax2, 'hot');
end
if n_fail > 0
    hold(ax2, 'on');
    scatter([results(fail_idx).cx_outer], [results(fail_idx).cy_outer], ...
            40, 'rx', 'LineWidth', 1.5);
    hold(ax2, 'off');
    legend(ax2, 'Pass (offset)', 'Fail (rejected)', ...
           'TextColor', 'w', 'Color', [0.15 0.15 0.2]);
end
set(ax2, 'Color', [0.15 0.15 0.2], 'XColor','w', 'YColor','w');
title(ax2, 'Ring offset per star (px) - green = used', 'Color', 'w');
xlabel(ax2, 'X (px)', 'Color', 'w'); ylabel(ax2, 'Y (px)', 'Color', 'w');
axis(ax2, 'image'); set(ax2, 'YDir', 'reverse');

ax3 = subplot(2,3,6);
if n_pass > 1
    histogram(ax3, offsets_pix_pass, max(5, round(n_pass/2)), ...
              'FaceColor', [0.3 0.6 1], 'EdgeColor', 'none');
elseif n_pass == 1
    bar(ax3, offsets_pix_pass, 'FaceColor', [0.3 0.6 1]);
end
xline(ax3, OFFSET_WARN_PIX, '--', 'Color', [1 0.7 0], 'LineWidth', 1.5, ...
      'Label', 'Warn');
xline(ax3, OFFSET_FAIL_PIX, '--', 'Color', [1 0.3 0.1], 'LineWidth', 1.5, ...
      'Label', 'Fail');
set(ax3, 'Color', [0.15 0.15 0.2], 'XColor', 'w', 'YColor', 'w');
xlabel(ax3, 'Ring center offset (px)', 'Color', 'w');
ylabel(ax3, 'Count', 'Color', 'w');
title(ax3, sprintf('Offset distribution (%d/%d stars used)', n_pass, nstars), ...
      'Color', 'w');

annotation('textbox', [0.02 0.01 0.96 0.05], ...
           'String', ['Verdict: ' verdict ...
                      sprintf('  |  Stars used: %d / %d', n_pass, nstars)], ...
           'Color', v_color, 'FontSize', 12, 'FontWeight', 'bold', ...
           'EdgeColor', 'none', 'BackgroundColor', 'none', ...
           'HorizontalAlignment', 'center');

% Save figure 1
fig1_base = fullfile(RESULT_DIR, 'collimation_analysis');
savefig(fig1, [fig1_base '.fig']);
exportgraphics(fig1, [fig1_base '.png'], 'Resolution', 150);
fprintf('Saved: %s.fig / .png\n', fig1_base);

% --- Figure 2: Annulus Quality Metrics
fig2 = figure('Name', 'Annulus Quality Metrics', 'NumberTitle', 'off', ...
       'Color', [0.12 0.12 0.15], 'Position', [50 620 1200 320]);

qc_fields  = {'contrast_ratio', 'annular_snr', 'hole_fill_frac', ...
              'ring_completeness', 'blob_circularity'};
qc_labels  = {'Contrast ratio', 'Annular SNR', 'Hole fill frac', ...
              'Ring completeness', 'Blob circularity'};
qc_thresh  = [MIN_CONTRAST_RATIO, MIN_ANNULAR_SNR, MAX_HOLE_FILL_FRAC, ...
              MIN_RING_COMPLETENESS, MIN_BLOB_CIRCULARITY];
qc_invert  = [false, false, true, false, false];

nq = numel(qc_fields);
for qi = 1:nq
    ax_q = subplot(1, nq, qi);
    vals = arrayfun(@(r) r.qc.(qc_fields{qi}), results);
    vals_disp = min(vals, 3 * qc_thresh(qi));

    bar_colors = zeros(nstars, 3);
    for k = 1:nstars
        if results(k).qc.pass
            bar_colors(k,:) = [0.3 0.6 1.0];
        else
            if qc_invert(qi)
                culprit = vals(k) > qc_thresh(qi);
            else
                culprit = vals(k) < qc_thresh(qi);
            end
            if culprit
                bar_colors(k,:) = [1.0 0.3 0.2];
            else
                bar_colors(k,:) = [0.6 0.6 0.6];
            end
        end
    end
    b = bar(ax_q, 1:nstars, vals_disp, 'FaceColor', 'flat');
    b.CData = bar_colors;
    yline(ax_q, qc_thresh(qi), '--y', 'LineWidth', 1.5);
    set(ax_q, 'Color', [0.15 0.15 0.2], 'XColor','w', 'YColor','w');
    title(ax_q, qc_labels{qi}, 'Color', 'w', 'FontSize', 9);
    xlabel(ax_q, 'Star #', 'Color', 'w');
    if qi == 1, ylabel(ax_q, 'Value', 'Color', 'w'); end
end

fig2_base = fullfile(RESULT_DIR, 'quality_metrics');
savefig(fig2, [fig2_base '.fig']);
exportgraphics(fig2, [fig2_base '.png'], 'Resolution', 150);
fprintf('Saved: %s.fig / .png\n', fig2_base);

% --- Figure 3: Polar Misalignment Vectors
fig3 = figure('Name', 'Misalignment Vectors', 'NumberTitle', 'off', ...
       'Color', [0.12 0.12 0.15], 'Position', [1470 50 500 500]);
ax4 = polaraxes;
set(ax4, 'Color', [0.15 0.15 0.2], ...
         'ThetaColor', 'w', 'RColor', 'w', 'GridColor', [0.5 0.5 0.5]);
hold(ax4, 'on');

cmap_pass = jet(max(n_pass, 1));
for qi = 1:n_pass
    k = pass_idx(qi);
    polarplot(ax4, deg2rad(results(k).angle_offset), results(k).offset_pix, ...
              'o', 'Color', cmap_pass(qi,:), 'MarkerSize', 8, ...
              'MarkerFaceColor', cmap_pass(qi,:));
end
for k = fail_idx
    polarplot(ax4, deg2rad(results(k).angle_offset), results(k).offset_pix, ...
              'd', 'Color', [0.5 0.5 0.5], 'MarkerSize', 5);
end
polarplot(ax4, [0, deg2rad(mean_angle)], [0, mean_offset], '-w', 'LineWidth', 2.5);
polarplot(ax4, deg2rad(mean_angle), mean_offset, ...
          'p', 'Color', 'w', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
title(ax4, sprintf('Misalignment vectors  (%d used, %d rejected)', ...
                   n_pass, n_fail), 'Color', 'w', 'FontSize', 11);
hold(ax4, 'off');

fig3_base = fullfile(RESULT_DIR, 'misalignment_vectors');
savefig(fig3, [fig3_base '.fig']);
exportgraphics(fig3, [fig3_base '.png'], 'Resolution', 150);
fprintf('Saved: %s.fig / .png\n', fig3_base);

fprintf('\nAnalysis complete. Adjust the secondary mirror opposite to\n');
fprintf('the mean misalignment direction (%.1f deg in image coords).\n', ...
        mean_angle);
fprintf('Quality filter kept %d of %d stars.\n', n_pass, nstars);
fprintf('\nAll results written to: %s\n\n', RESULT_DIR);

%% =========================================================================
%  LOCAL FUNCTIONS
%  =========================================================================

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