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
%  OUTPUTS
%    - Annotated image: detected stars with ellipse overlays
%    - Per-star table:  centroid offset, ellipticity, radial symmetry score
%    - Collimation verdict: GOOD / WARNING / POOR with directional hint
%    Written by Claud Sonnet 4.6 April 2026

%

clear all;
clc;
close all;

%% =========================================================================
%  PARAMETERS  –  edit these for your telescope and camera
%  =========================================================================
% Leave '' to open a file-selection dialog
FITS_PATH         = '/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST_20260412';
FITS_FILE         = fullfile(FITS_PATH,'LAST.01.01.01_20260412.173436.632_clear_GoodFocusTest_001_001_001_sci_raw_Image_1.fits');       

% Detector & optics
PIXEL_SCALE       = 1.25;      % arcsec / pixel (for your plate scale)
OBSTRUCTION_RATIO = 0.40;     % Secondary mirror diameter / Primary diameter
                              % (e.g. 0.35 for a 35% linear obstruction)

% Star detection
BG_SIGMA_CLIP     = 3.0;      % Sigma-clipping iterations for background
DETECT_THRESH     = 300.0;      % Detection threshold in sigma above background
STAR_RADIUS_MIN   = 8;        % Minimum expected star radius (pixels)
STAR_RADIUS_MAX   = 50;       % Maximum expected star radius (pixels)
MAX_STARS         = 100;       % Analyse at most this many brightest stars

% Ellipse fitting
MIN_RING_PIXELS   = 60;       % Minimum pixels on each ring for a valid fit

% Collimation thresholds
OFFSET_WARN_PIX   = 1.5;      % Ring center offset (px) triggering WARNING
OFFSET_FAIL_PIX   = 4.0;      % Ring center offset (px) triggering POOR
ELLIP_WARN        = 0.05;     % Ellipticity excess (1 - b/a) for WARNING
ELLIP_FAIL        = 0.12;     % Ellipticity excess for POOR

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
raw   = double(fitsread(FITS_FILE));          % Primary image HDU

% Handle multi-extension files or 3-D cubes (take first plane / HDU)
if ndims(raw) == 3
    raw = raw(:,:,1);
    warning('3-D array detected – using first plane only.');
end

[nrows, ncols] = size(raw);
fprintf('Image size: %d x %d pixels\n', ncols, nrows);

%% =========================================================================
%  2. PREPROCESSING
%  =========================================================================
% --- 2a. Debayer / grayscale (if needed) ---
% If your sensor is monochrome this section does nothing.
% For colour FITS just convert to luminance:
if size(raw, 3) == 3
    raw = 0.2126*raw(:,:,1) + 0.7152*raw(:,:,2) + 0.0722*raw(:,:,3);
end

% --- 2b. Sigma-clipped sky background and noise estimation ---
img_flat = raw(:);
bg_med   = median(img_flat);
bg_std   = std(img_flat);
for k = 1:5                                   % iterative sigma clip
    mask     = abs(img_flat - bg_med) < BG_SIGMA_CLIP * bg_std;
    bg_med   = median(img_flat(mask));
    bg_std   = std(img_flat(mask));
end
fprintf('Sky background: %.1f  Noise (sigma): %.2f\n', bg_med, bg_std);

% --- 2c. Subtract background ---
img = raw - bg_med;

%% =========================================================================
%  3. STAR DETECTION
%  =========================================================================
% Annular (donut) stars have a dark hole in the centre, so we cannot use
% a simple centroid finder. Instead we:
%   (a) threshold the image to find bright ring blobs
%   (b) use morphological operations to fill the central hole
%   (c) measure connected-region properties to find star candidates

thresh_map = img > DETECT_THRESH * bg_std;

% Fill holes so each donut becomes a filled disk
thresh_fill = imfill(thresh_map, 'holes');

% Remove small noise blobs
min_area = pi * STAR_RADIUS_MIN^2;
max_area = pi * STAR_RADIUS_MAX^2;
thresh_clean = bwareaopen(thresh_fill, round(min_area));

% Label connected components
CC   = bwconncomp(thresh_clean);
rp   = regionprops(CC, img, ...
       'Centroid','Area','MajorAxisLength','MinorAxisLength', ...
       'Orientation','WeightedCentroid','BoundingBox','Solidity');

% Filter by area range
valid = ([rp.Area] > min_area) & ([rp.Area] < max_area);
rp    = rp(valid);

if isempty(rp)
    error(['No stars detected. Try adjusting DETECT_THRESH, ' ...
           'STAR_RADIUS_MIN/MAX, or check the FITS file.']);
end

% Sort by brightness (sum of pixel values in bounding box region) descending
flux = zeros(1, numel(rp));
for k = 1:numel(rp)
    bb  = round(rp(k).BoundingBox);
    r1  = max(1, bb(2)); r2 = min(nrows, bb(2)+bb(4));
    c1  = max(1, bb(1)); c2 = min(ncols, bb(1)+bb(3));
    flux(k) = sum(sum(img(r1:r2, c1:c2)));
end
[~, sidx] = sort(flux, 'descend');
rp   = rp(sidx(1:min(MAX_STARS, numel(sidx))));
nstars = numel(rp);
fprintf('Stars selected for analysis: %d\n', nstars);

%% =========================================================================
%  4. RING EXTRACTION AND ELLIPSE FITTING
%  =========================================================================
% For each detected star we:
%   (i)  Cut a sub-image (stamp) around the star
%   (ii) Detect the outer ring edge (bright ring) and inner ring edge (hole)
%   (iii) Fit ellipses to both edges using least-squares conic fitting

% results = struct('id',{}, 'cx_outer',{}, 'cy_outer',{}, ...
%                  'cx_inner',{}, 'cy_inner',{}, ...
%                  'offset_pix',{}, 'offset_arcsec',{}, ...
%                  'ellipticity_outer',{}, 'ellipticity_inner',{}, ...
%                  'radial_variance',{}, 'angle_offset',{});

% results = struct([]);   % empty — fields are defined on first assignment
for k = 1:nstars
    bb  = rp(k).BoundingBox;
    pad = round(0.3 * max(bb(3), bb(4)));           % 30% padding
    r1  = max(1, round(bb(2)) - pad);
    r2  = min(nrows, round(bb(2)+bb(4)) + pad);
    c1  = max(1, round(bb(1)) - pad);
    c2  = min(ncols, round(bb(1)+bb(3)) + pad);

    stamp    = img(r1:r2, c1:c2);                   % sub-image (stamp)
    [sr, sc] = size(stamp);

    % -- Outer ring: Canny edge of thresholded stamp
    stamp_norm = mat2gray(stamp);
    stamp_filled = imfill(stamp_norm,'holes');
    edges_out  = edge(stamp_filled, 'canny', [0.05 0.15]);

    % -- Inner ring: detect the boundary of the central dark hole
    % Threshold inverted stamp to find the dark interior
    hole = stamp_filled-stamp_norm;
    hole_mask  = hole > (0.15 * max(hole(:)));
    % hole_mask1  = imfill(hole_mask0, 'holes');
    % hole_mask2  = bwareaopen(hole_mask1, 20);
    edges_in   = bwperim(hole_mask);

    % Fit ellipses via algebraic (least-squares) conic method
    [ry_o, rx_o] = find(edges_out);
    [ry_i, rx_i] = find(edges_in);

    res_k = struct();
    res_k.id = k;

    % Outer ellipse
    if numel(rx_o) >= MIN_RING_PIXELS
        ep_o = fit_ellipse_algebraic(rx_o, ry_o);
        res_k.cx_outer         = ep_o.cx + c1 - 1;
        res_k.cy_outer         = ep_o.cy + r1 - 1;
        res_k.ellipticity_outer = 1 - ep_o.minor/ep_o.major;
        res_k.ep_outer_stamp   = ep_o;
    else
        % Fall back to region-prop centroid
        ep_o = [];
        res_k.cx_outer = rp(k).WeightedCentroid(1);
        res_k.cy_outer = rp(k).WeightedCentroid(2);
        res_k.ellipticity_outer = 1 - rp(k).MinorAxisLength / ...
                                       max(rp(k).MajorAxisLength, 1);
        res_k.ep_outer_stamp = [];
        warning('Star %d: insufficient outer ring pixels for ellipse fit.', k);
    end

    % Inner ellipse
    if numel(rx_i) >= MIN_RING_PIXELS
        ep_i = fit_ellipse_algebraic(rx_i, ry_i);
        res_k.cx_inner         = ep_i.cx + c1 - 1;
        res_k.cy_inner         = ep_i.cy + r1 - 1;
        res_k.ellipticity_inner = 1 - ep_i.minor/ep_i.major;
        res_k.ep_inner_stamp   = ep_i;
    else
        ep_i = [];
        res_k.cx_inner = rp(k).WeightedCentroid(1);
        res_k.cy_inner = rp(k).WeightedCentroid(2);
        res_k.ellipticity_inner = res_k.ellipticity_outer;
        res_k.ep_inner_stamp = [];
        warning('Star %d: insufficient inner ring pixels for ellipse fit.', k);
    end

    % -- Ring center offset (key collimation metric)
    dx = res_k.cx_outer - res_k.cx_inner;
    dy = res_k.cy_outer - res_k.cy_inner;
    res_k.offset_pix    = sqrt(dx^2 + dy^2);
    res_k.offset_arcsec = res_k.offset_pix * PIXEL_SCALE;
    res_k.angle_offset  = atan2d(dy, dx);   % direction of misalignment

    % -- Radial brightness variance (azimuthal symmetry check)
    cx_s = res_k.cx_outer - c1 + 1;
    cy_s = res_k.cy_outer - r1 + 1;
    res_k.radial_variance = radial_brightness_variance(stamp, cx_s, cy_s, ...
                                   round(ep_o.minor * 0.6), ...
                                   round(ep_o.major * 1.1));

    results(k) = res_k;
end

%% =========================================================================
%  5. COLLIMATION ASSESSMENT
%  =========================================================================
offsets_pix   = [results.offset_pix];
ellipticities = [results.ellipticity_outer];

mean_offset    = mean(offsets_pix);
median_offset  = median(offsets_pix);
mean_ellip     = mean(ellipticities);

% Mean misalignment vector (tells you WHICH WAY to adjust)
mean_dx = mean([results.cx_outer] - [results.cx_inner]);
mean_dy = mean([results.cy_outer] - [results.cy_inner]);
mean_angle = atan2d(mean_dy, mean_dx);

% Verdict
if mean_offset < OFFSET_WARN_PIX && mean_ellip < ELLIP_WARN
    verdict = 'GOOD   – collimation is within acceptable limits';
    v_color = [0.0 0.6 0.2];
elseif mean_offset < OFFSET_FAIL_PIX && mean_ellip < ELLIP_FAIL
    verdict = 'WARNING – slight collimation error detected';
    v_color = [0.8 0.5 0.0];
else
    verdict = 'POOR   – significant collimation error; adjustment recommended';
    v_color = [0.8 0.1 0.1];
end

%% =========================================================================
%  6. RESULTS TABLE
%  =========================================================================
fprintf('\n==============================================================\n');
fprintf('  COLLIMATION REPORT\n');
fprintf('==============================================================\n');
fprintf('  File  : %s\n', FITS_FILE);
fprintf('  Stars : %d analysed\n', nstars);
fprintf('  Mean ring-centre offset : %.2f px  (%.2f arcsec)\n', ...
         mean_offset, mean_offset * PIXEL_SCALE);
fprintf('  Median offset           : %.2f px\n', median_offset);
fprintf('  Mean outer ellipticity  : %.3f  (0 = perfect circle)\n', mean_ellip);
fprintf('  Mean misalignment angle : %.1f deg (image coords)\n', mean_angle);
fprintf('\n  VERDICT: %s\n', verdict);
fprintf('==============================================================\n\n');

fprintf('  %-5s %-10s %-12s %-14s %-10s %-12s\n', ...
        'Star', 'Offset(px)', 'Offset(arcsec)', 'Ellipticity', ...
        'RadVar', 'Angle(deg)');
fprintf('  %s\n', repmat('-', 1, 66));
for k = 1:nstars
    fprintf('  %-5d %-10.2f %-12.2f %-14.3f %-10.4f %-12.1f\n', ...
            k, results(k).offset_pix, results(k).offset_arcsec, ...
            results(k).ellipticity_outer, results(k).radial_variance, ...
            results(k).angle_offset);
end
fprintf('\n');

%% =========================================================================
%  7. DIAGNOSTIC PLOTS
%  =========================================================================
figure('Name', 'Collimation Analysis', 'NumberTitle', 'off', ...
       'Color', [0.12 0.12 0.15], 'Position', [50 50 1400 900]);

% --- 7a. Full annotated field ---
ax1 = subplot(2,3,[1,2,4,5]);
display_img = log1p(max(img,0));               % log stretch for display
imshow(display_img, [], 'Parent', ax1);
colormap(ax1, 'gray');
hold(ax1, 'on');
title(ax1, 'Annotated star field', 'Color', 'w', 'FontSize', 11);

cmap_heat = jet(nstars);
for k = 1:nstars
    clr = cmap_heat(k,:);
    % Outer ellipse
    if ~isempty(results(k).ep_outer_stamp)
        draw_ellipse(ax1, results(k).cx_outer, results(k).cy_outer, ...
                     results(k).ep_outer_stamp.major, ...
                     results(k).ep_outer_stamp.minor, ...
                     results(k).ep_outer_stamp.angle, clr, 1.5);
    end
    % Inner ellipse
    if ~isempty(results(k).ep_inner_stamp)
        draw_ellipse(ax1, results(k).cx_inner, results(k).cy_inner, ...
                     results(k).ep_inner_stamp.major, ...
                     results(k).ep_inner_stamp.minor, ...
                     results(k).ep_inner_stamp.angle, clr, 1.0);
    end
    % Offset vector
    plot(ax1, [results(k).cx_outer, results(k).cx_inner], ...
              [results(k).cy_outer, results(k).cy_inner], ...
              'Color', [1 0.4 0.1], 'LineWidth', 1.2);
    text(ax1, results(k).cx_outer+5, results(k).cy_outer-5, ...
         sprintf('%d', k), 'Color', clr, 'FontSize', 8);
end
hold(ax1, 'off');

% --- 7b. Ring center offsets (scatter) ---
ax2 = subplot(2,3,3);
scatter([results.cx_outer], [results.cy_outer], 60, offsets_pix, 'filled');
colorbar; colormap(ax2, 'hot');
set(ax2, 'Color', [0.15 0.15 0.2], 'XColor','w', 'YColor','w');
title(ax2, 'Ring offset per star (px)', 'Color', 'w');
xlabel(ax2, 'X (px)', 'Color', 'w'); ylabel(ax2, 'Y (px)', 'Color', 'w');
axis(ax2, 'image'); set(ax2, 'YDir', 'reverse');

% --- 7c. Offset histogram ---
ax3 = subplot(2,3,6);
histogram(ax3, offsets_pix, max(5, round(nstars/2)), ...
          'FaceColor', [0.3 0.6 1], 'EdgeColor', 'none');
xline(ax3, OFFSET_WARN_PIX, '--', 'Color', [1 0.7 0], 'LineWidth', 1.5, ...
      'Label', 'Warn');
xline(ax3, OFFSET_FAIL_PIX, '--', 'Color', [1 0.3 0.1], 'LineWidth', 1.5, ...
      'Label', 'Fail');
set(ax3, 'Color', [0.15 0.15 0.2], 'XColor', 'w', 'YColor', 'w');
xlabel(ax3, 'Ring center offset (px)', 'Color', 'w');
ylabel(ax3, 'Count', 'Color', 'w');
title(ax3, 'Offset distribution', 'Color', 'w');

% Overall verdict text
annotation('textbox', [0.02 0.01 0.96 0.05], ...
           'String', ['Verdict: ' verdict], ...
           'Color', v_color, 'FontSize', 12, 'FontWeight', 'bold', ...
           'EdgeColor', 'none', 'BackgroundColor', 'none', ...
           'HorizontalAlignment', 'center');

% --- 7d. Polar plot of misalignment vectors ---
figure('Name', 'Misalignment Vectors', 'NumberTitle', 'off', ...
       'Color', [0.12 0.12 0.15], 'Position', [1470 50 500 500]);
ax4 = polaraxes;
set(ax4, 'Color', [0.15 0.15 0.2], ...
         'ThetaColor', 'w', 'RColor', 'w', 'GridColor', [0.5 0.5 0.5]);
hold(ax4, 'on');
for k = 1:nstars
    polarplot(ax4, deg2rad(results(k).angle_offset), results(k).offset_pix, ...
              'o', 'Color', cmap_heat(k,:), 'MarkerSize', 8, ...
              'MarkerFaceColor', cmap_heat(k,:));
end
% Mean vector
polarplot(ax4, [0, deg2rad(mean_angle)], [0, mean_offset], ...
          '-w', 'LineWidth', 2.5);
polarplot(ax4, deg2rad(mean_angle), mean_offset, ...
          'p', 'Color', 'w', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
title(ax4, 'Misalignment vectors', 'Color', 'w', 'FontSize', 11);
hold(ax4, 'off');

fprintf('Analysis complete. Adjust the secondary mirror opposite to\n');
fprintf('the mean misalignment direction (%.1f deg in image coords).\n\n', ...
        mean_angle);

%% =========================================================================
%  LOCAL FUNCTIONS
%  =========================================================================

function ep = fit_ellipse_algebraic(x, y)
%FIT_ELLIPSE_ALGEBRAIC  Least-squares algebraic ellipse fit (Fitzgibbon 1999)
%  Returns struct with fields: cx, cy, major, minor, angle (degrees)
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
    cx = (C*D - B*F) / denom;
    cy = (A*F - B*D) / denom;
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
%  Low variance = symmetric (well-collimated)
%  High variance = asymmetric brightness (coma / tilt)
    [rows, cols] = size(stamp);
    [X, Y] = meshgrid(1:cols, 1:rows);
    R = sqrt((X - cx).^2 + (Y - cy).^2);
    ring = stamp(R >= r_in & R <= r_out);
    if isempty(ring)
        rv = 0; return;
    end
    rv = std(ring) / (mean(ring) + eps);   % coefficient of variation
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