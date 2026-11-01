function S = analyze_annulus_cutout(Stamp,Args)
% ANALYZE_ANNULUS_CUTOUT  Extract collimation statistics from a single
%   pre-cropped image cutout containing one annular (donut) star image.
%
%  The cutout may be larger than the annulus (the star occupies only part
%  of it).  Step 2 locates the annulus bounding box automatically by
%  thresholding, filling the central hole with imfill, and finding the
%  largest connected blob.  All geometry is derived from that tight
%  sub-stamp and then reported back in full-stamp pixel coordinates.
%
%  INPUT
%    Stamp   - 2-D numeric array (double or integer), raw or bg-subtracted.
%              Set Args.SubtractBackground = true for raw data.
%
%    Args    - (optional) name-value pairs or struct with any of the fields
%              listed below.
%
%  ARGS FIELDS  (all optional - defaults shown)
%    Args.PixelScale          = 1.25    % arcsec / pixel
%    Args.ObstructionRatio    = 0.40    % secondary / primary diameter ratio
%    Args.MinRingPixels       = 60      % min edge pixels for valid ellipse fit
%    Args.SubtractBackground  = false   % estimate & subtract bg before processing
%    Args.BgSigmaClip         = 3.0     % sigma-clip half-width for corner bg
%    Args.BgNoise             = []      % externally measured bg noise (sigma);
%                                       % if empty, estimated from stamp corners
%    Args.DetectNSigma        = 3.0     % threshold for annulus detection
%                                       % (multiples of bg_std above bg_med)
%    Args.MinBlobArea         = 100     % minimum filled-blob area (px^2) to
%                                       % accept as a real annulus
%    Args.SectorFillThreshold = 0.5     % fraction of ring_median a sector must
%                                       % exceed to count as "lit"
%    Args.NSectors            = 16      % azimuthal sectors for completeness check
%    Args.MinContrastRatio    = 2.0     % quality threshold
%    Args.MinAnnularSnr       = 6.0     % quality threshold
%    Args.MaxHoleFillFrac     = 0.45    % quality threshold
%    Args.MinRingCompleteness = 0.65    % quality threshold
%    Args.MinBlobCircularity  = 0.30    % quality threshold
%    Args.DiagnosticPlot      = true    % show 6-panel diagnostic figure;
%                                       % set false to suppress
%
%  OUTPUT  S - struct with fields:
%
%   Bounding-box (location of annulus inside the input Stamp, 1-based)
%    S.BB_r1, S.BB_r2   - row range of the tight annulus sub-stamp
%    S.BB_c1, S.BB_c2   - col range of the tight annulus sub-stamp
%    S.SubStamp         - the tight sub-image that was analysed
%
%   Background / noise  (estimated from full Stamp corners)
%    S.BgMedian         - background level (ADU)
%    S.BgStd            - background noise (ADU)
%
%   Geometry  (all coordinates in FULL STAMP pixels, 1-based)
%    S.CxOuter          - outer ellipse centre x
%    S.CyOuter          - outer ellipse centre y
%    S.MajorOuter       - outer ellipse semi-major axis (px)
%    S.MinorOuter       - outer ellipse semi-minor axis (px)
%    S.AngleOuter       - outer ellipse position angle (deg, from x-axis)
%    S.OuterRadius      - mean outer radius = (MajorOuter+MinorOuter)/2
%    S.CxInner          - inner ellipse centre x
%    S.CyInner          - inner ellipse centre y
%    S.MajorInner       - inner ellipse semi-major axis (px)
%    S.MinorInner       - inner ellipse semi-minor axis (px)
%    S.AngleInner       - inner ellipse position angle (deg)
%    S.HoleArea         - area of the central dark hole (px^2), from imfill
%    S.RingArea         - area of the bright ring pixels (px^2)
%    S.FilledArea       - total filled-blob area  = HoleArea + RingArea (px^2)
%
%   Collimation metrics
%    S.Dx               - signed x-component of outer-inner offset (px)
%    S.Dy               - signed y-component
%    S.OffsetPix        - |offset|  (px)
%    S.OffsetNorm       - OffsetPix / OuterRadius  (dimensionless)
%    S.OffsetArcsec     - offset in arcseconds
%    S.AngleOffset      - direction of inner->outer shift (deg, image coords)
%
%   Photometric / symmetry
%    S.EllipticityOuter - 1 - MinorOuter/MajorOuter  (0 = circle)
%    S.EllipticityInner - 1 - MinorInner/MajorInner
%    S.RadialVariance   - azimuthal CoV of ring brightness (low = symmetric)
%
%   Annulus quality metrics
%    S.ContrastRatio    - ring_median / hole_median
%    S.AnnularSnr       - ring_mean / bg_std
%    S.HoleFillFrac     - hole_median / ring_median
%    S.RingCompleteness - fraction of azimuthal sectors that are illuminated
%    S.BlobCircularity  - filled_blob_area / bounding_circle_area
%
%   Quality gate
%    S.QualityPass      - logical: true if all five metrics pass
%    S.FailReasons      - cell array of strings naming failed metrics
%
%   Raw ellipse-fit structs (sub-stamp coords; empty when fallback used)
%    S.EpOuter, S.EpInner
%
%  Written by Claude Sonnet 4.6, June 2026

arguments
    Stamp                              (:,:) {mustBeNumeric, mustBeReal}
    Args.PixelScale          (1,1) double  = 1.25
    Args.ObstructionRatio    (1,1) double  = 0.40
    Args.MinRingPixels       (1,1) double  = 60
    Args.SubtractBackground  (1,1) logical = false
    Args.BgSigmaClip         (1,1) double  = 3.0
    Args.BgNoise             (:,1) double  = []
    Args.DetectNSigma        (1,1) double  = 3.0
    Args.MinBlobArea         (1,1) double  = 100
    Args.SectorFillThreshold (1,1) double  = 0.5
    Args.NSectors            (1,1) double  = 16
    Args.MinContrastRatio    (1,1) double  = 2.0
    Args.MinAnnularSnr       (1,1) double  = 6.0
    Args.MaxHoleFillFrac     (1,1) double  = 0.45
    Args.MinRingCompleteness (1,1) double  = 0.65
    Args.MinBlobCircularity  (1,1) double  = 0.30
    Args.DiagnosticPlot      (1,1) logical = true   % set false to suppress figure
end

Stamp = double(Stamp);
[sr, sc] = size(Stamp);
% smmothing
% Stamp = imdiffusefilt(Stamp);
% Stamp = wiener2(Stamp, [5 5]);
% Stamp = imgaussfilt(Stamp,2);
Stamp = medfilt2(Stamp,[2,2]);

% =========================================================================
%  STEP 1 – Background estimation from stamp corners
% =========================================================================
dr = max(1, round(0.10 * sr));
dc = max(1, round(0.10 * sc));
corner_pixels = [ Stamp(1:dr,       1:dc      );
                  Stamp(1:dr,       sc-dc+1:sc);
                  Stamp(sr-dr+1:sr, 1:dc      );
                  Stamp(sr-dr+1:sr, sc-dc+1:sc) ];
corner_pixels = corner_pixels(:);

bg_med = median(corner_pixels);
bg_std = std(corner_pixels);
for iter = 1:5
    mask_bg = abs(corner_pixels - bg_med) < Args.BgSigmaClip * bg_std;
    if sum(mask_bg) < 4, break; end
    bg_med = median(corner_pixels(mask_bg));
    bg_std = std(corner_pixels(mask_bg));
end
bg_std = max(bg_std, 1.0);

if ~isempty(Args.BgNoise)
    bg_std = Args.BgNoise;
end

if Args.SubtractBackground
    Stamp        = Stamp - bg_med;
    bg_med_report = 0;
else
    bg_med_report = bg_med;
end

% =========================================================================
%  STEP 2 – Locate the annulus inside the cutout
%
%  Strategy:
%    a) Threshold the background-subtracted stamp at DetectNSigma * bg_std
%       to get a binary mask of the bright ring pixels.
%    b) Fill the central hole with imfill so the donut becomes a solid disk.
%    c) Keep only blobs larger than MinBlobArea; pick the largest one.
%    d) Derive:
%         - tight bounding box  -> SubStamp (the annulus sub-image)
%         - filled-blob area    -> FilledArea
%         - hole mask           -> HoleArea  (filled minus original ring)
%         - ring area           -> RingArea
%    e) All downstream geometry uses SubStamp coordinates, then offsets
%       are added back to report positions in full-stamp coordinates.
% =========================================================================

% --- 2a. Threshold (work on bg-subtracted values) -----------------------
thresh_level = Args.DetectNSigma * bg_std;
if Args.SubtractBackground
    ring_binary = Stamp > thresh_level;
else
    ring_binary = (Stamp - bg_med) > thresh_level;
end

% --- 2b. Fill holes to turn donut into solid disk -----------------------
filled_binary = imfill(ring_binary, 'holes');

% --- 2c. Pick the largest filled blob -----------------------------------
filled_binary = bwareaopen(filled_binary, round(Args.MinBlobArea));
CC   = bwconncomp(filled_binary);

if CC.NumObjects == 0
    error(['analyze_annulus_cutout: no annulus blob found above threshold. ' ...
          'Try lowering Args.DetectNSigma (currently %.1f).'], Args.DetectNSigma);
end

rp_blob = regionprops(CC, 'Area', 'BoundingBox');
[~, biggest] = max([rp_blob.Area]);
bb_blob      = rp_blob(biggest).BoundingBox;   % [x0 y0 width height] (regionprops convention)

% --- 2d. Areas -----------------------------------------------------------
% Rebuild mask for just the chosen blob
chosen_mask   = false(sr, sc);
chosen_mask(CC.PixelIdxList{biggest}) = true;

ring_mask     = chosen_mask & ring_binary;     % bright ring pixels only
hole_mask_bin = chosen_mask & ~ring_binary;    % dark interior pixels

FilledArea = sum(chosen_mask(:));
RingArea   = sum(ring_mask(:));
HoleArea   = sum(hole_mask_bin(:));

% --- 2e. Tight bounding box -> SubStamp ---------------------------------
% regionprops BoundingBox: [x0, y0, width, height] where x0,y0 are the
% top-left corner - 0.5 (sub-pixel).  Convert to integer row/col ranges.
bb_c1 = max(1,  floor(bb_blob(1)));
bb_r1 = max(1,  floor(bb_blob(2)));
bb_c2 = min(sc, ceil( bb_blob(1) + bb_blob(3)));
bb_r2 = min(sr, ceil( bb_blob(2) + bb_blob(4)));

SubStamp      = Stamp(bb_r1:bb_r2, bb_c1:bb_c2);
sub_chosen    = chosen_mask(bb_r1:bb_r2, bb_c1:bb_c2);
sub_ring      = ring_mask(bb_r1:bb_r2, bb_c1:bb_c2);
sub_hole      = hole_mask_bin(bb_r1:bb_r2, bb_c1:bb_c2);

[ss_r, ss_c] = size(SubStamp);

% Column/row offsets to convert SubStamp coords -> full-Stamp coords
col_offset = bb_c1 - 1;   % cx_full = cx_sub + col_offset
row_offset = bb_r1 - 1;   % cy_full = cy_sub + row_offset

% Fallback radius estimates from the filled blob (used if ellipse fit fails)
r_outer_est = 0.5 * max(bb_blob(3), bb_blob(4));
r_inner_est = r_outer_est * Args.ObstructionRatio;

% =========================================================================
%  STEP 3 – Annulus quality metrics  (computed on SubStamp)
% =========================================================================
[Xg, Yg] = meshgrid(1:ss_c, 1:ss_r);

% Centroid of the filled blob in SubStamp coords (area-weighted)
sub_cx = mean(Xg(sub_chosen));
sub_cy = mean(Yg(sub_chosen));

ring_pixels  = SubStamp(sub_ring);
hole_pixels  = SubStamp(sub_hole);

ring_median  = median(ring_pixels);
ring_mean    = mean(ring_pixels);
hole_median  = median(hole_pixels);  % sub_hole may be empty if hole is tiny

% 1. Contrast ratio
if isempty(hole_pixels) || hole_median <= 0
    contrast_ratio = Inf;
else
    contrast_ratio = ring_median / hole_median;
end

% 2. Annular SNR
annular_snr = ring_mean / (bg_std + eps);

% 3. Hole fill fraction
if ring_median > 0
    hole_fill_frac = hole_median / ring_median;
else
    hole_fill_frac = 1.0;
end

% 4. Ring completeness (azimuthal sector coverage around blob centroid)
n_sectors    = Args.NSectors;
sector_angle = 2*pi / n_sectors;
theta_map    = atan2(Yg - sub_cy, Xg - sub_cx);
sector_lit   = false(1, n_sectors);
for s = 1:n_sectors
    a0 = -pi + (s-1) * sector_angle;
    a1 = a0 + sector_angle;
    sec_mask = sub_ring & (theta_map >= a0) & (theta_map < a1);
    if any(sec_mask(:)) && ...
       median(SubStamp(sec_mask)) > Args.SectorFillThreshold * ring_median
        sector_lit(s) = true;
    end
end
ring_completeness = sum(sector_lit) / n_sectors;

% 5. Blob circularity  (filled area / bounding-circle area)
bounding_r         = 0.5 * sqrt(bb_blob(3)^2 + bb_blob(4)^2);
bounding_circ_area = pi * bounding_r^2;
blob_circularity   = FilledArea / bounding_circ_area;
%==========================================================================
%******************* Here is the crucial code
% =========================================================================
%  STEP 4 – Ellipse fitting on outer and inner ring edges
%           (everything in SubStamp coordinates)
% =========================================================================
% Normalise SubStamp for edge detection (mat2gray handles the range)
sub_norm   = mat2gray(SubStamp);
% Re-fill to ensure the outer boundary is fully closed before Canny
sub_filled = imfill(sub_norm, 'holes');

% Outer edge: Canny on the filled normalised sub-stamp
edges_out = edge(sub_filled, 'canny', [0.05 0.15]);

% Inner edge: boundary of the central dark hole
%   hole_diff highlights where imfill added pixels (i.e. the interior)
hole_diff  = sub_filled - sub_norm;
hole_mask2 = hole_diff > (0.15 * max(hole_diff(:)));
edges_in   = bwperim(hole_mask2);

[ry_o, rx_o] = find(edges_out);
[ry_i, rx_i] = find(edges_in);

% --- Outer ellipse (SubStamp coords) ------------------------------------
if numel(rx_o) >= Args.MinRingPixels
    ep_o = fit_ellipse_algebraic(rx_o, ry_o);
    cx_outer_sub  = ep_o.cx;
    cy_outer_sub  = ep_o.cy;
    major_outer   = ep_o.major;
    minor_outer   = ep_o.minor;
    angle_outer   = ep_o.angle;
    ellip_outer   = 1 - ep_o.minor / max(ep_o.major, eps);
    outer_radius  = (ep_o.major + ep_o.minor) / 2;
    ep_outer_out  = ep_o;
else
    ep_o = struct('cx', sub_cx, 'cy', sub_cy, ...
                  'major', r_outer_est, 'minor', r_inner_est, 'angle', 0);
    cx_outer_sub  = sub_cx;
    cy_outer_sub  = sub_cy;
    major_outer   = r_outer_est;
    minor_outer   = r_inner_est;
    angle_outer   = 0;
    ellip_outer   = 1 - r_inner_est / max(r_outer_est, eps);
    outer_radius  = (r_outer_est + r_inner_est) / 2;
    ep_outer_out  = [];
    warning('analyze_annulus_cutout: outer edge only %d px (<%d); fallback used.', ...
            numel(rx_o), Args.MinRingPixels);
end

% --- Inner ellipse (SubStamp coords) ------------------------------------
if numel(rx_i) >= Args.MinRingPixels
    ep_i = fit_ellipse_algebraic(rx_i, ry_i);
    cx_inner_sub  = ep_i.cx;
    cy_inner_sub  = ep_i.cy;
    major_inner   = ep_i.major;
    minor_inner   = ep_i.minor;
    angle_inner   = ep_i.angle;
    ellip_inner   = 1 - ep_i.minor / max(ep_i.major, eps);
    ep_inner_out  = ep_i;
else
    cx_inner_sub  = cx_outer_sub;
    cy_inner_sub  = cy_outer_sub;
    major_inner   = r_inner_est;
    minor_inner   = r_inner_est;
    angle_inner   = 0;
    ellip_inner   = ellip_outer;
    ep_inner_out  = [];
    warning('analyze_annulus_cutout: inner edge only %d px (<%d); fallback used.', ...
            numel(rx_i), Args.MinRingPixels);
end

% --- Map centres back to full-stamp coordinates -------------------------
cx_outer = cx_outer_sub + col_offset;
cy_outer = cy_outer_sub + row_offset;
cx_inner = cx_inner_sub + col_offset;
cy_inner = cy_inner_sub + row_offset;

% =========================================================================
%  STEP 5 – Collimation metrics
% =========================================================================
dx            = cx_outer - cx_inner;
dy            = cy_outer - cy_inner;
offset_pix    = sqrt(dx^2 + dy^2);
offset_norm   = offset_pix / max(outer_radius, eps);
offset_arcsec = offset_pix * Args.PixelScale;
angle_offset  = atan2d(dy, dx);

% =========================================================================
%  STEP 6 – Radial brightness variance (azimuthal symmetry, on SubStamp)
% =========================================================================
radial_var = radial_brightness_variance(SubStamp, cx_outer_sub, cy_outer_sub, ...
                 round(ep_o.minor * 0.6), ...
                 round(ep_o.major * 1.1));

% =========================================================================
%  STEP 7 – Quality gate
% =========================================================================
fail_contrast    = contrast_ratio    < Args.MinContrastRatio;
fail_snr         = annular_snr       < Args.MinAnnularSnr;
fail_hole        = hole_fill_frac    > Args.MaxHoleFillFrac;
fail_complete    = ring_completeness < Args.MinRingCompleteness;
fail_circularity = blob_circularity  < Args.MinBlobCircularity;

quality_pass = ~(fail_contrast || fail_snr || fail_hole || ...
                 fail_complete || fail_circularity);

fail_reasons = {};
if fail_contrast,    fail_reasons{end+1} = 'contrast';     end
if fail_snr,         fail_reasons{end+1} = 'SNR';          end
if fail_hole,        fail_reasons{end+1} = 'holefill';     end
if fail_complete,    fail_reasons{end+1} = 'completeness'; end
if fail_circularity, fail_reasons{end+1} = 'circularity';  end

% =========================================================================
%  STEP 8 – Pack output struct
% =========================================================================
% Bounding box
S.BB_r1        = bb_r1;
S.BB_r2        = bb_r2;
S.BB_c1        = bb_c1;
S.BB_c2        = bb_c2;
S.SubStamp     = SubStamp;

% Background
S.BgMedian     = bg_med_report;
S.BgStd        = bg_std;

% Outer ellipse (full-stamp coords)
S.CxOuter      = cx_outer;
S.CyOuter      = cy_outer;
S.MajorOuter   = major_outer;
S.MinorOuter   = minor_outer;
S.AngleOuter   = angle_outer;
S.OuterRadius  = outer_radius;

% Inner ellipse (full-stamp coords)
S.CxInner      = cx_inner;
S.CyInner      = cy_inner;
S.MajorInner   = major_inner;
S.MinorInner   = minor_inner;
S.AngleInner   = angle_inner;

% Areas derived from imfill
S.FilledArea   = FilledArea;
S.RingArea     = RingArea;
S.HoleArea     = HoleArea;

% Collimation
S.Dx           = dx;
S.Dy           = dy;
S.OffsetPix    = offset_pix;
S.OffsetNorm   = offset_norm;
S.OffsetArcsec = offset_arcsec;
S.AngleOffset  = angle_offset;

% Photometric
S.EllipticityOuter = ellip_outer;
S.EllipticityInner = ellip_inner;
S.RadialVariance   = radial_var;

% Quality
S.ContrastRatio    = contrast_ratio;
S.AnnularSnr       = annular_snr;
S.HoleFillFrac     = hole_fill_frac;
S.RingCompleteness = ring_completeness;
S.BlobCircularity  = blob_circularity;
S.QualityPass      = quality_pass;
S.FailReasons      = fail_reasons;

% Raw ellipse structs (SubStamp coords; empty if fallback used)
S.EpOuter      = ep_outer_out;
S.EpInner      = ep_inner_out;

% =========================================================================
%  STEP 9 – Diagnostic figure  (silenced when Args.DiagnosticPlot = false)
% =========================================================================
if Args.DiagnosticPlot
    BG_FIG = [0.12 0.12 0.15];
    BG_AX  = [0.15 0.15 0.20];

    fig = figure('Name', 'analyze_annulus_cutout  diagnostic', ...
                 'NumberTitle', 'off', 'Color', BG_FIG, ...
                 'Position', [60 60 1300 780]);

    % Shared log-stretch of the full stamp for display
    disp_full = log1p(max(Stamp, 0));
    disp_sub  = log1p(max(SubStamp, 0));

    % ------------------------------------------------------------------
    %  Panel 1 – Full cutout (log stretch) with bounding box
    % ------------------------------------------------------------------
    ax1 = subplot(2, 3, 1);
    imshow(disp_full, [], 'Parent', ax1);
    colormap(ax1, 'gray');
    hold(ax1, 'on');
    % Draw bounding box as a rectangle
    rectangle('Parent', ax1, ...
              'Position', [bb_c1, bb_r1, bb_c2-bb_c1, bb_r2-bb_r1], ...
              'EdgeColor', [1 0.8 0], 'LineWidth', 1.8, 'LineStyle', '--');
    title(ax1, 'Full cutout + annulus bounding box', ...
          'Color', 'w', 'FontSize', 9);
    hold(ax1, 'off');
    set(ax1, 'Color', BG_AX);

    % ------------------------------------------------------------------
    %  Panel 2 – Threshold binary (ring pixels before hole-fill)
    % ------------------------------------------------------------------
    ax2 = subplot(2, 3, 2);
    % ring_binary is in full-stamp coords; show the sub-region
    imshow(ring_binary(bb_r1:bb_r2, bb_c1:bb_c2), 'Parent', ax2);
    title(ax2, sprintf('Ring binary  (%.1f sigma threshold)', Args.DetectNSigma), ...
          'Color', 'w', 'FontSize', 9);
    set(ax2, 'Color', BG_AX);

    % ------------------------------------------------------------------
    %  Panel 3 – Filled blob with hole mask overlaid
    % ------------------------------------------------------------------
    ax3 = subplot(2, 3, 3);
    % Build an RGB overlay: filled blob = grey, hole = red, ring = white
    sub_filled_show = chosen_mask(bb_r1:bb_r2, bb_c1:bb_c2);
    sub_hole_show   = hole_mask_bin(bb_r1:bb_r2, bb_c1:bb_c2);
    sub_ring_show   = ring_mask(bb_r1:bb_r2, bb_c1:bb_c2);
    rgb = zeros(ss_r, ss_c, 3);
    rgb(:,:,1) = 0.55 * sub_filled_show + 0.45 * sub_ring_show;  % R
    rgb(:,:,2) = 0.55 * sub_filled_show;                          % G
    rgb(:,:,3) = 0.55 * sub_filled_show;                          % B
    % Colour the hole red
    rgb(:,:,1) = rgb(:,:,1) + 0.7 * sub_hole_show;
    rgb(:,:,2) = rgb(:,:,2) - 0.3 * double(sub_hole_show);
    rgb(:,:,3) = rgb(:,:,3) - 0.3 * double(sub_hole_show);
    rgb = max(0, min(1, rgb));
    imshow(rgb, 'Parent', ax3);
    title(ax3, sprintf('Filled blob (grey) + hole (red)  |  FilledArea=%d  HoleArea=%d', ...
                       FilledArea, HoleArea), ...
          'Color', 'w', 'FontSize', 9);
    set(ax3, 'Color', BG_AX);

    % ------------------------------------------------------------------
    %  Panel 4 – SubStamp with detected edge points
    % ------------------------------------------------------------------
    ax4 = subplot(2, 3, 4);
    imshow(disp_sub, [], 'Parent', ax4);
    colormap(ax4, 'gray');
    hold(ax4, 'on');
    % Outer edge points (cyan dots)
    plot(ax4, rx_o, ry_o, '.', 'Color', [0.2 0.9 1.0], ...
         'MarkerSize', 2);
    % Inner edge points (orange dots)
    plot(ax4, rx_i, ry_i, '.', 'Color', [1.0 0.5 0.1], ...
         'MarkerSize', 2);
    legend(ax4, sprintf('Outer (%d px)', numel(rx_o)), ...
                sprintf('Inner (%d px)', numel(rx_i)), ...
           'TextColor', 'w', 'Color', BG_AX, 'EdgeColor', 'none', ...
           'FontSize', 7, 'Location', 'southoutside');
    title(ax4, 'SubStamp + edge points', 'Color', 'w', 'FontSize', 9);
    hold(ax4, 'off');
    set(ax4, 'Color', BG_AX);

    % ------------------------------------------------------------------
    %  Panel 5 – SubStamp with fitted ellipses and offset vector
    % ------------------------------------------------------------------
    ax5 = subplot(2, 3, 5);
    imshow(disp_sub, [], 'Parent', ax5);
    colormap(ax5, 'gray');
    hold(ax5, 'on');

    % Outer ellipse (cyan)
    draw_ellipse_diag(ax5, cx_outer_sub, cy_outer_sub, ...
                      major_outer, minor_outer, angle_outer, ...
                      [0.2 0.9 1.0], 1.8);
    plot(ax5, cx_outer_sub, cy_outer_sub, '+', ...
         'Color', [0.2 0.9 1.0], 'MarkerSize', 8, 'LineWidth', 1.5);

    % Inner ellipse (orange)
    draw_ellipse_diag(ax5, cx_inner_sub, cy_inner_sub, ...
                      major_inner, minor_inner, angle_inner, ...
                      [1.0 0.5 0.1], 1.4);
    plot(ax5, cx_inner_sub, cy_inner_sub, '+', ...
         'Color', [1.0 0.5 0.1], 'MarkerSize', 8, 'LineWidth', 1.5);

    % Offset vector (yellow arrow, outer->inner)
    if offset_pix > 0.5
        quiver(ax5, cx_outer_sub, cy_outer_sub, ...
               cx_inner_sub - cx_outer_sub, cy_inner_sub - cy_outer_sub, ...
               0, 'Color', [1 0.95 0.1], 'LineWidth', 2.0, ...
               'MaxHeadSize', 1.5);
    end

    title(ax5, sprintf('Fitted ellipses  |  offset=%.2fpx (%.4f norm)  angle=%.1fdeg', ...
                       offset_pix, offset_norm, angle_offset), ...
          'Color', 'w', 'FontSize', 9);
    hold(ax5, 'off');
    set(ax5, 'Color', BG_AX);

    % ------------------------------------------------------------------
    %  Panel 6 – Quality metrics text summary
    % ------------------------------------------------------------------
    ax6 = subplot(2, 3, 6);
    set(ax6, 'Color', BG_AX, 'XColor', 'none', 'YColor', 'none');
    axis(ax6, 'off');

    % Build a table of metric | value | threshold | pass/fail
    metrics = { ...
        'Contrast ratio',    contrast_ratio,    Args.MinContrastRatio,    ~fail_contrast,    'min'; ...
        'Annular SNR',       annular_snr,        Args.MinAnnularSnr,       ~fail_snr,         'min'; ...
        'Hole fill frac',    hole_fill_frac,     Args.MaxHoleFillFrac,     ~fail_hole,        'max'; ...
        'Ring completeness', ring_completeness,  Args.MinRingCompleteness, ~fail_complete,    'min'; ...
        'Blob circularity',  blob_circularity,   Args.MinBlobCircularity,  ~fail_circularity, 'min'; ...
    };

    % Column header
    text(ax6, 0.02, 0.97, ...
         sprintf('%-20s  %7s  %9s  %s', 'Metric', 'Value', 'Threshold', 'Status'), ...
         'Units', 'normalized', 'Color', [0.85 0.85 0.85], ...
         'FontName', 'Monospaced', 'FontSize', 8.5, 'FontWeight', 'bold', ...
         'VerticalAlignment', 'top');
    text(ax6, 0.02, 0.90, repmat('-', 1, 48), ...
         'Units', 'normalized', 'Color', [0.5 0.5 0.5], ...
         'FontName', 'Monospaced', 'FontSize', 8.5, ...
         'VerticalAlignment', 'top');

    row_y = 0.83;
    for mi = 1:size(metrics, 1)
        name_  = metrics{mi, 1};
        val_   = metrics{mi, 2};
        thr_   = metrics{mi, 3};
        pass_  = metrics{mi, 4};
        sense_ = metrics{mi, 5};  % 'min' or 'max'
        if pass_
            status_str = 'PASS';
            clr        = [0.3 0.9 0.3];
        else
            status_str = 'FAIL';
            clr        = [1.0 0.3 0.2];
        end
        if isinf(val_)
            val_str = sprintf('%7s', 'Inf');
        else
            val_str = sprintf('%7.3f', val_);
        end
        line_str = sprintf('%-20s  %s  %5s %5.3f  %s', ...
                           name_, val_str, sense_, thr_, status_str);
        text(ax6, 0.02, row_y, line_str, ...
             'Units', 'normalized', 'Color', clr, ...
             'FontName', 'Monospaced', 'FontSize', 8.5, ...
             'VerticalAlignment', 'top');
        row_y = row_y - 0.10;
    end

    % Overall verdict
    row_y = row_y - 0.04;
    text(ax6, 0.02, row_y, repmat('-', 1, 48), ...
         'Units', 'normalized', 'Color', [0.5 0.5 0.5], ...
         'FontName', 'Monospaced', 'FontSize', 8.5, ...
         'VerticalAlignment', 'top');
    row_y = row_y - 0.09;
    if quality_pass
        vstr = 'PASS';  vclr = [0.3 0.9 0.3];
    else
        reasons_str = strjoin(fail_reasons, ', ');
        vstr = sprintf('FAIL  (%s)', reasons_str);
        vclr = [1.0 0.3 0.2];
    end
    text(ax6, 0.02, row_y, sprintf('Overall: %s', vstr), ...
         'Units', 'normalized', 'Color', vclr, ...
         'FontName', 'Monospaced', 'FontSize', 9, 'FontWeight', 'bold', ...
         'VerticalAlignment', 'top');

    % Extra stats below
    row_y = row_y - 0.13;
    extra = { ...
        sprintf('Outer radius : %.1f px', outer_radius); ...
        sprintf('Ellip outer  : %.4f', ellip_outer); ...
        sprintf('Ellip inner  : %.4f', ellip_inner); ...
        sprintf('Radial var   : %.4f', radial_var); ...
        sprintf('Angle offset : %.1f deg', angle_offset); ...
        sprintf('HoleArea     : %d px^2', HoleArea); ...
        sprintf('RingArea     : %d px^2', RingArea); ...
    };
    for ei = 1:numel(extra)
        text(ax6, 0.02, row_y, extra{ei}, ...
             'Units', 'normalized', 'Color', [0.75 0.85 1.0], ...
             'FontName', 'Monospaced', 'FontSize', 8.5, ...
             'VerticalAlignment', 'top');
        row_y = row_y - 0.09;
    end

    title(ax6, 'Quality metrics', 'Color', 'w', 'FontSize', 9);

    sgtitle(fig, ...
            sprintf('analyze\\_annulus\\_cutout  |  offset=%.2fpx  norm=%.4f  angle=%.1fdeg', ...
                    offset_pix, offset_norm, angle_offset), ...
            'Color', 'w', 'FontSize', 11);
end  % if DiagnosticPlot

end  % function analyze_annulus_cutout


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function ep = fit_ellipse_algebraic(x, y)
    x = x(:); y = y(:);
    D = [x.^2, x.*y, y.^2, x, y, ones(size(x))];
    S = D' * D;
    C = zeros(6,6); C(1,3) = 2; C(3,1) = 2; C(2,2) = -1;
    [V, lam] = eig(S \ C);
    lambda   = diag(lam);
    pos_def  = (lambda > 1e-10) & isfinite(lambda);
    if ~any(pos_def), ep = fallback_ellipse(x, y); return; end
    [~, idx] = min(lambda(pos_def));
    vp = V(:, pos_def);
    a  = vp(:, idx);
    ep = conic_to_geometric(a);
end

function ep = conic_to_geometric(a)
    A = a(1); B = a(2)/2; C = a(3); D = a(4)/2; F = a(5)/2; G = a(6);
    denom = B^2 - A*C;
    if abs(denom) < 1e-12, ep = fallback_ellipse([], []); return; end
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
    if isempty(x)
        ep.cx = 0; ep.cy = 0; ep.major = 10; ep.minor = 10; ep.angle = 0;
    else
        ep.cx    = mean(x);  ep.cy    = mean(y);
        ep.major = (max(x) - min(x)) / 2;
        ep.minor = (max(y) - min(y)) / 2;
        ep.angle = 0;
    end
end

function rv = radial_brightness_variance(stamp, cx, cy, r_in, r_out)
    [rows, cols] = size(stamp);
    [X, Y] = meshgrid(1:cols, 1:rows);
    R    = sqrt((X - cx).^2 + (Y - cy).^2);
    ring = stamp(R >= r_in & R <= r_out);
    if isempty(ring), rv = 0; return; end
    rv = std(ring) / (mean(ring) + eps);
end

function draw_ellipse_diag(ax, cx, cy, major, minor, angle_deg, color, lw)
%DRAW_ELLIPSE_DIAG  Overlay a single ellipse on axes ax.
    t   = linspace(0, 2*pi, 300);
    ang = deg2rad(angle_deg);
    xe  = major * cos(t);
    ye  = minor * sin(t);
    xr  = cx + xe * cos(ang) - ye * sin(ang);
    yr  = cy + xe * sin(ang) + ye * cos(ang);
    plot(ax, xr, yr, 'Color', color, 'LineWidth', lw);
end