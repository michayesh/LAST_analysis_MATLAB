% OutOfFocusCollimationAnalysis
% 1. Read the fits image file
% 2. Preprocess it using pipeline.generic.prePrep
% 3. Debias flatten the image using pipeline.generic.basicCalib and generic
%    dark and flat images
% 4. Divide the image into 24 non-overlapping subImages (crops)
% 5. Loop over the crops, use imUtil.sources.findRings to detect annuli,
%    pick the brightest, cut it out, and run analyze_annulus_cutout on it.
% 6. Collect per-subimage results into SubResults struct array.
% 7. Compute aggregate collimation statistics across all subimages and
%    report verdict + save plots.
%
% Micha, June 2026 + Claude Sonnet 4.6

clear;
close all;

%% =========================================================================
%  PARAMETERS
%% =========================================================================
THRESHOLD       = 6.;
RINGS_TEMP_BANK = [11 20; 15 33; 21 50; 31 60; 50 100];
MAXRAD          = 100;      % half-size of annulus cutout (pixels)
EDGE_GUARD      = MAXRAD;   % reject rings within this many px of sub-image edge

% Collimation verdict thresholds (normalised offset = offset / outer_radius)
OFFSET_WARN_NORM = 0.05;
OFFSET_FAIL_NORM = 0.15;
ELLIP_WARN       = 0.05;
ELLIP_FAIL       = 0.12;

% Pixel scale for the instrument
PIXEL_SCALE = 1.25;   % arcsec / pixel

FITS_PATH  = '/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST.01.06.01_20260517/';
ImageName  = 'LAST.01.06.01_20260517.183147.208_clear_OOF1000_010_001_001_sci_raw_Image_1.fits';
dark_file  = ['/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST.01.04.04_20260430/calib/' ...
              'LAST.01.04.03_20241223.230027.415_clear_____dark_proc_Image_1.fits'];
flat_file  = ['/home/micha/Dropbox/WAO/LAST_analysis/Out_of_focus_Images/LAST.01.04.04_20260430/calib/' ...
              'LAST.01.04.03_20241204.150643.935_clear_____twflat_proc_Image_1.fits'];

%% =========================================================================
%  1. READ IMAGE
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
%  2. PREPROCESS
%% =========================================================================
Args.prePrepArgs = {};
% Do not check PSF 
AI = pipeline.generic.prePrep(ai, Args.prePrepArgs{:}, ...
                               'GlobalBackLevel', false, 'MinNim', 1, 'GlobalBadPSF', false);

%% =========================================================================
%  3. CALIBRATE
%% =========================================================================
Args.basicCalibArgs = {};
AI = pipeline.generic.basicCalib(AI, Args.basicCalibArgs{:}, CI, 'UpdateJD', false);

%% =========================================================================
%  4. DIVIDE INTO SUB-IMAGES
%% =========================================================================
croppedImage = AI.Image(1:9576, 1:6388);
ImageSize    = [6388 9576];
SubSize      = [1597 1596];
[SubImageCCDSEC, ~, ~, ~, ~] = imUtil.cut.gridSubImage(ImageSize, SubSize);
nSubIm       = size(SubImageCCDSEC, 1);

%% =========================================================================
%  5 + 6.  LOOP: DETECT RINGS, PICK BRIGHTEST, ANALYSE CUTOUT
%% =========================================================================
% Pre-allocate SubResults as a struct array.
% Fields added per iteration:
%   .SubIndex         - subimage index (1..nSubIm)
%   .CCDSEC           - [r1 r2 c1 c2] of this subimage in the full image
%   .NRingsFound      - total rings detected before edge/threshold cuts
%   .NRingsAboveMedian- rings above 50th-percentile amplitude threshold
%   .BrightestPosSubIm- [x y] of brightest ring in subimage coords
%   .BrightestPosFullIm-[x y] of brightest ring in full-image coords
%   .BrightestAmplitude- findRings amplitude of the brightest ring
%   .Cutout           - 2-D double array: the annulus stamp
%   .Stats            - struct returned by analyze_annulus_cutout
%   .Valid            - logical: false if no usable ring was found

SubResults = repmat(struct( ...
    'SubIndex',           [], ...
    'CCDSEC',             [], ...
    'NRingsFound',        0,  ...
    'NRingsAboveMedian',  0,  ...
    'BrightestPosSubIm',  [], ...
    'BrightestPosFullIm', [], ...
    'BrightestAmplitude', NaN, ...
    'Cutout',             [], ...
    'Stats',              [], ...
    'Valid',              false), nSubIm, 1);

% Args passed to analyze_annulus_cutout (uses same pixel scale / obstruction)
% cutoutArgs.PixelScale       = PIXEL_SCALE;
% cutoutArgs.ObstructionRatio = 0.40;

% for isub = 1:nSubIm
for isub = [10 11 14 15] % only the central group

    sec       = SubImageCCDSEC(isub, :);   % [r1 r2 c1 c2]
    SubImage  = imUtil.cut.trim2d(croppedImage, sec);
    [SizeY, SizeX] = size(SubImage);

    SubResults(isub).SubIndex = isub;
    SubResults(isub).CCDSEC   = sec;

    % --- Detect rings -------------------------------------------------------
    [Pos, ~, ~] = imUtil.sources.findRings(SubImage, ...
                    'CalcBack',    true, ...
                    'Threshold',   THRESHOLD, ...
                    'FunBackVar',  @imUtil.background.modeVar_SampleHist, ...
                    'RingsTempBank', RINGS_TEMP_BANK);

    % Pos columns: [x, y, amplitude, ...]
    if isempty(Pos)
        fprintf('SubImage %2d: no rings detected.\n', isub);
        continue
    end

    % --- Remove rings too close to the sub-image edge -----------------------
    edgeFlag = Pos(:,1) > EDGE_GUARD          & ...
               Pos(:,2) > EDGE_GUARD          & ...
               Pos(:,1) < (SizeX - EDGE_GUARD) & ...
               Pos(:,2) < (SizeY - EDGE_GUARD);
    Pos = Pos(edgeFlag, :);

    if isempty(Pos)
        fprintf('SubImage %2d: all rings rejected by edge guard.\n', isub);
        continue
    end

    SubResults(isub).NRingsFound = size(Pos, 1);

    % Rings above 50th-percentile amplitude (stored for future use)
    ampThresh = prctile(Pos(:, 3), 50);
    SubResults(isub).NRingsAboveMedian = sum(Pos(:, 3) >= ampThresh);

    fprintf('SubImage %2d: %d rings found, %d above median amplitude.\n', ...
            isub, SubResults(isub).NRingsFound, SubResults(isub).NRingsAboveMedian);

    % --- Pick the brightest ring --------------------------------------------
    [maxAmp, maxIdx]  = max(Pos(:, 3));
    xBright           = Pos(maxIdx, 1);   % x in subimage coords (col)
    yBright           = Pos(maxIdx, 2);   % y in subimage coords (row)

    SubResults(isub).BrightestAmplitude = maxAmp;
    SubResults(isub).BrightestPosSubIm  = [xBright, yBright];

    % --- Convert to full-image coordinates ----------------------------------
    % sec = [r1 r2 c1 c2]: subimage rows are sec(1):sec(2), cols are sec(3):sec(4)
    % findRings x = column direction, y = row direction
    xFull = xBright + sec(3) - 1;   % column offset
    yFull = yBright + sec(1) - 1;   % row offset
    SubResults(isub).BrightestPosFullIm = [xFull, yFull];

    % --- Extract annulus cutout (MAXRAD half-size box) ----------------------
    % Cutout row range: yBright +/- MAXRAD  (clamp to subimage bounds)
    rMin = max(1,     floor(yBright - MAXRAD));
    rMax = min(SizeY, ceil( yBright + MAXRAD));
    cMin = max(1,     floor(xBright - MAXRAD));
    cMax = min(SizeX, ceil( xBright + MAXRAD));

    annulusCutout = double(SubImage(rMin:rMax, cMin:cMax));
    SubResults(isub).Cutout = annulusCutout;

    % --- Analyse the cutout -------------------------------------------------
    try
        % S = analyze_annulus_cutout(annulusCutout, cutoutArgs);
        S = analyze_annulus_cutout(annulusCutout,'SubtractBackground',true,'DetectNSigma',6);
        SubResults(isub).Stats = S;
        SubResults(isub).Valid = true;
    catch ME
        warning('SubImage %d: analyze_annulus_cutout failed: %s', isub, ME.message);
    end

end   % for isub

%% =========================================================================
%  7.  AGGREGATE COLLIMATION STATISTICS  (valid subimages only)
%% =========================================================================
validMask  = [SubResults.Valid];
validIdx   = find(validMask);
nValid     = numel(validIdx);

fprintf('\n--- Aggregate results: %d / %d subimages valid ---\n', nValid, nSubIm);

if nValid == 0
    warning('No valid subimages – cannot compute collimation statistics.');
else
    % Extract per-subimage scalar metrics from Stats structs
    offsets_pix   = arrayfun(@(r) r.Stats.OffsetPix,        SubResults(validIdx));
    offsets_norm  = arrayfun(@(r) r.Stats.OffsetNorm,       SubResults(validIdx));
    ellipticities = arrayfun(@(r) r.Stats.EllipticityOuter, SubResults(validIdx));
    dx_vals       = arrayfun(@(r) r.Stats.Dx,               SubResults(validIdx));
    dy_vals       = arrayfun(@(r) r.Stats.Dy,               SubResults(validIdx));
    angles        = arrayfun(@(r) r.Stats.AngleOffset,      SubResults(validIdx));
    outer_radii   = arrayfun(@(r) r.Stats.OuterRadius,      SubResults(validIdx));
    radial_vars   = arrayfun(@(r) r.Stats.RadialVariance,   SubResults(validIdx));
    qual_pass     = arrayfun(@(r) r.Stats.QualityPass,      SubResults(validIdx));

    % Global positions of analysed stars (for spatial offset map)
    pos_full      = vertcat(SubResults(validIdx).BrightestPosFullIm);  % [nValid x 2]

    % Summary statistics
    mean_outer_rad    = mean(outer_radii);
    std_outer_rad     = std(outer_radii);
    mean_offset_pix   = mean(offsets_pix);
    std_offset_pix    = std(offsets_pix);
    median_offset_pix = median(offsets_pix);
    mean_offset_norm  = mean(offsets_norm);
    std_offset_norm   = std(offsets_norm);
    mean_ellip        = mean(ellipticities);
    std_ellip         = std(ellipticities);
    mean_dx           = mean(dx_vals);
    mean_dy           = mean(dy_vals);
    std_dx            = std(dx_vals);
    std_dy            = std(dy_vals);
    mean_angle        = atan2d(mean_dy, mean_dx);
    std_angle         = std(angles);

    % Verdict
    if mean_offset_norm < OFFSET_WARN_NORM && mean_ellip < ELLIP_WARN
        verdict = 'GOOD - collimation is within acceptable limits';
        v_color = [0.0 0.6 0.2];
    elseif mean_offset_norm < OFFSET_FAIL_NORM && mean_ellip < ELLIP_FAIL
        verdict = 'WARNING - slight collimation error detected';
        v_color = [0.8 0.5 0.0];
    else
        verdict = 'POOR - significant collimation error; adjustment recommended';
        v_color = [0.8 0.1 0.1];
    end

    % -----------------------------------------------------------------------
    %  PRINT REPORT
    % -----------------------------------------------------------------------
    sep = repmat('=', 1, 66);
    fprintf('\n%s\n', sep);
    fprintf('  COLLIMATION REPORT\n');
    fprintf('%s\n', sep);
    fprintf('  File            : %s\n', ImageName);
    fprintf('  SubImages       : %d total, %d valid, %d passed quality filter\n', ...
            nSubIm, nValid, sum(qual_pass));
    fprintf('  Mean outer radius       : %.2f +/- %.2f px\n', mean_outer_rad, std_outer_rad);
    fprintf('  Mean offset             : %.2f +/- %.2f px  (%.3f +/- %.3f norm)\n', ...
            mean_offset_pix, std_offset_pix, mean_offset_norm, std_offset_norm);
    fprintf('  Median offset           : %.2f px\n', median_offset_pix);
    fprintf('  Mean offset (arcsec)    : %.3f +/- %.3f\n', ...
            mean_offset_pix * PIXEL_SCALE, std_offset_pix * PIXEL_SCALE);
    fprintf('  Mean ellipticity        : %.3f +/- %.3f\n', mean_ellip, std_ellip);
    fprintf('  Mean misalignment angle : %.1f +/- %.1f deg\n', mean_angle, std_angle);
    fprintf('\n  VERDICT: %s\n', verdict);
    fprintf('%s\n\n', sep);

    % Per-subimage table
    fprintf('  %-6s %-6s %-12s %-10s %-12s %-10s %-10s\n', ...
            'SubIm', 'Valid', 'Offset(px)', 'OffNorm', 'Ellipticity', 'Angle', 'QualPass');
    fprintf('  %s\n', repmat('-', 1, 68));
    for qi = 1:nValid
        k = validIdx(qi);
        fprintf('  %-6d %-6s %-12.2f %-10.4f %-12.3f %-10.1f %-10s\n', ...
                k, 'yes', offsets_pix(qi), offsets_norm(qi), ...
                ellipticities(qi), angles(qi), ...
                mat2str(qual_pass(qi)));
    end
    % Pad failed subimages at the bottom
    for k = find(~validMask)
        fprintf('  %-6d %-6s\n', k, 'no');
    end
    fprintf('\n');

    % -----------------------------------------------------------------------
    %  PLOTS
    % -----------------------------------------------------------------------
    BG_FIG = [0.12 0.12 0.15];
    BG_AX  = [0.15 0.15 0.20];

    % -- Fig 1: Spatial offset map over full detector -----------------------
    fig1 = figure('Name', 'Ring Offset Map', 'Color', BG_FIG, ...
                  'Position', [50 50 800 600]);
    ax1  = axes(fig1);
    set(ax1, 'Color', BG_AX, 'XColor', 'w', 'YColor', 'w');
    hold(ax1, 'on');

    % Draw subimage grid for reference
    for isub = 1:nSubIm
        sec = SubImageCCDSEC(isub, :);
        % sec = [r1 r2 c1 c2]  ->  plot box in (x=col, y=row) space
        rx = [sec(3) sec(4) sec(4) sec(3) sec(3)];
        ry = [sec(1) sec(1) sec(2) sec(2) sec(1)];
        plot(ax1, rx, ry, 'Color', [0.35 0.35 0.4], 'LineWidth', 0.5);
    end

    sc_h = scatter(ax1, pos_full(:,1), pos_full(:,2), ...
                   60, offsets_pix, 'filled');
    cb = colorbar(ax1);  cb.Color = 'w';
    cb.Label.String = 'Ring centre offset (px)';  cb.Label.Color = 'w';
    colormap(ax1, 'hot');

    % Overlay misalignment vectors scaled for visibility
    scale = 500;   % pixels per unit offset_pix for arrow length
    for qi = 1:nValid
        quiver(ax1, pos_full(qi,1), pos_full(qi,2), ...
               dx_vals(qi) * scale, dy_vals(qi) * scale, ...
               0, 'Color', [0.3 0.8 1], 'LineWidth', 1.2, ...
               'MaxHeadSize', 0.5);
    end

    title(ax1, sprintf('Offset map  |  %d/%d subimages  |  Verdict: %s', ...
                       nValid, nSubIm, verdict), ...
          'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');
    xlabel(ax1, 'X full image (px)', 'Color', 'w');
    ylabel(ax1, 'Y full image (px)', 'Color', 'w');
    set(ax1, 'YDir', 'reverse');
    hold(ax1, 'off');

    % -- Fig 2: Offset histogram (normalised) --------------------------------
    fig2 = figure('Name', 'Offset Histogram', 'Color', BG_FIG, ...
                  'Position', [200 100 700 450]);
    ax2  = axes(fig2);
    set(ax2, 'Color', BG_AX, 'XColor', 'w', 'YColor', 'w');
    hold(ax2, 'on');
    if nValid > 1
        histogram(ax2, offsets_norm, max(5, round(nValid/2)), ...
                  'FaceColor', [0.2 0.8 0.6], 'EdgeColor', 'none');
    else
        bar(ax2, offsets_norm, 'FaceColor', [0.2 0.8 0.6]);
    end
    xline(ax2, OFFSET_WARN_NORM, '--', 'Color', [1 0.7 0], 'LineWidth', 1.8, ...
          'Label', sprintf('Warn (%.3f)', OFFSET_WARN_NORM), 'LabelColor', [1 0.7 0]);
    xline(ax2, OFFSET_FAIL_NORM, '--', 'Color', [1 0.3 0.1], 'LineWidth', 1.8, ...
          'Label', sprintf('Fail (%.3f)', OFFSET_FAIL_NORM), 'LabelColor', [1 0.3 0.1]);
    xline(ax2, mean_offset_norm, '-', 'Color', [0.4 1 0.4], 'LineWidth', 2.0, ...
          'Label', sprintf('Mean (%.4f)', mean_offset_norm), 'LabelColor', [0.4 1 0.4]);
    xlabel(ax2, 'Normalised offset  (offset / outer\_radius)', 'Color', 'w');
    ylabel(ax2, 'Count', 'Color', 'w');
    title(ax2, sprintf('Offset distribution (normalised)  |  %d subimages', nValid), ...
          'Color', 'w', 'FontSize', 11);
    hold(ax2, 'off');

    % -- Fig 3: Polar misalignment plot -------------------------------------
    fig3 = figure('Name', 'Misalignment Vectors', 'Color', BG_FIG, ...
                  'Position', [350 150 560 520]);
    ax3  = polaraxes('Parent', fig3);
    set(ax3, 'Color', BG_AX, 'ThetaColor', 'w', 'RColor', 'w', ...
             'GridColor', [0.5 0.5 0.5]);
    hold(ax3, 'on');
    ax3.RLim = [0 1];
    theta_c  = linspace(0, 2*pi, 361);
    polarplot(ax3, theta_c, repmat(OFFSET_WARN_NORM, 1, 361), ...
              '--', 'Color', [1 0.7 0], 'LineWidth', 1.2);
    polarplot(ax3, theta_c, repmat(OFFSET_FAIL_NORM, 1, 361), ...
              '--', 'Color', [1 0.3 0.1], 'LineWidth', 1.2);
    cmap_v = jet(max(nValid, 1));
    for qi = 1:nValid
        polarplot(ax3, deg2rad(angles(qi)), offsets_norm(qi), ...
                  'o', 'Color', cmap_v(qi,:), 'MarkerSize', 8, ...
                  'MarkerFaceColor', cmap_v(qi,:));
    end
    polarplot(ax3, [0, deg2rad(mean_angle)], [0, mean_offset_norm], ...
              '-w', 'LineWidth', 2.5);
    polarplot(ax3, deg2rad(mean_angle), mean_offset_norm, ...
              'p', 'Color', 'w', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    title(ax3, sprintf('Misalignment vectors  |  %d subimages', nValid), ...
          'Color', 'w', 'FontSize', 11);
    hold(ax3, 'off');

    % -- Fig 4: Cutout gallery (up to 24 stamps in a grid) ------------------
    nCols  = 6;
    nRows  = ceil(nValid / nCols);
    fig4   = figure('Name', 'Annulus Cutout Gallery', 'Color', BG_FIG, ...
                    'Position', [50 50 1300 220*nRows]);
    for qi = 1:nValid
        k      = validIdx(qi);
        ax_g   = subplot(nRows, nCols, qi);
        stamp  = SubResults(k).Cutout;
        imshow(log1p(max(stamp, 0)), [], 'Parent', ax_g);
        colormap(ax_g, 'gray');

        S = SubResults(k).Stats;
        title(ax_g, sprintf('Sub %d  off=%.2fpx  ang=%.0f', ...
              k, S.OffsetPix, S.AngleOffset), ...
              'Color', 'w', 'FontSize', 7);

        % Overlay outer and inner ellipses
        hold(ax_g, 'on');
        if ~isempty(S.EpOuter)
            draw_ellipse_local(ax_g, S.CxOuter, S.CyOuter, ...
                               S.MajorOuter, S.MinorOuter, S.AngleOuter, ...
                               [0.3 0.8 1], 1.2);
        end
        if ~isempty(S.EpInner)
            draw_ellipse_local(ax_g, S.CxInner, S.CyInner, ...
                               S.MajorInner, S.MinorInner, S.AngleInner, ...
                               [1 0.5 0.1], 1.0);
        end
        % Offset vector
        plot(ax_g, [S.CxOuter S.CxInner], [S.CyOuter S.CyInner], ...
             'Color', [1 0.9 0], 'LineWidth', 1.0);
        hold(ax_g, 'off');
        set(ax_g, 'XColor', 'none', 'YColor', 'none');
    end
    sgtitle(fig4, sprintf('Annulus cutouts  |  Verdict: %s', verdict), ...
            'Color', 'w', 'FontSize', 11, 'Interpreter', 'none');

    fprintf('Adjust secondary mirror OPPOSITE to mean misalignment direction: %.1f deg.\n\n', mean_angle);
end   % if nValid > 0


%% =========================================================================
%  LOCAL FUNCTION
%% =========================================================================
function draw_ellipse_local(ax, cx, cy, major, minor, angle_deg, color, lw)
    t   = linspace(0, 2*pi, 200);
    ang = deg2rad(angle_deg);
    xe  = major * cos(t);
    ye  = minor * sin(t);
    xr  = cx + xe * cos(ang) - ye * sin(ang);
    yr  = cy + xe * sin(ang) + ye * cos(ang);
    plot(ax, xr, yr, 'Color', color, 'LineWidth', lw);
end