%% RASA11_aberration_analysis.m
% Field-dependent aberration diagnosis of a RASA 11 (f/2.2) star-field image
% using AstroPack (https://github.com/EranOfek/AstroPack).
%
% What it measures, per isolated unsaturated star, then across the field:
%   * Coma        : (a) centroid shift  = broad-window centroid - core centroid
%                   (b) normalized spin-1 3rd moment  c = (M30+M12, M21+M03)/s^3
%                   Direction convention: vector points from the sharp core
%                   toward the flare. Radial component > 0 => flare points
%                   AWAY from the optical axis.
%   * Trefoil     : spin-3 3rd moment (M30-3M12, 3M21-M03)/s^3 (pinched optics)
%   * Astigmatism : spin-2 ellipticity (e1,e2) from window-deconvolved 2nd moments
%   * Focus       : FWHM over the field -> fitted best-focus surface
%                   (defocus offset, sensor tilt, field curvature)
%   * ePSF mosaic : oversampled effective PSF in a 3x3 grid of the field
%
% Field models fitted (field coords x,y in mm from the optical center):
%   coma  [Cx;Cy] = C0 + A*[x;y]          -> constant (decenter/collimation)
%                                             + linear (3rd-order field coma)
%                                             -> coma-free point = -A\C0
%   astig e1 = e10 + a1 x + a2 y + kA (x^2-y^2)/Rmax^2
%         e2 = e20 + b1 x + b2 y + kA (2xy)/Rmax^2
%   focus FWHM^2 = F0^2 + ((z0 + tx x + ty y + cz r^2)/N)^2      [um]
%
% Memory: the image is kept in its native class and processed crop by crop
% (P.CropGrid, with an overlap margin P.CropMargin). Each crop gets its own
% background and source detection; only per-star measurements and running
% ePSF sums are kept, so peak memory is set by one crop, not the full frame.
%
% Requirements: MATLAB R2020a+ (exportgraphics), AstroPack on the path.
% No toolboxes needed (robust fits are implemented locally).
%
% Coordinates follow the FITS/ds9 convention: X = column, Y = row, plots use
% "axis xy" so Y increases upward.

clear; close all;

%% ========================= User parameters =============================
P.ImageDir      = '/home/micha/Dropbox/WAO/LAST_analysis/psf_map/LAST.01.01.02';
P.ImageFile     = '';      % '' -> file dialog in P.ImageDir; else file name (in P.ImageDir) or full path
P.HDU           = 1;
P.OutDirPrefix  = 'IQA_';   % output folder = <prefix><image name>, created in the image's folder (IQA = Image Quality Analysis)
P.Satur         = [];      % ADU. [] -> header SATURVAL/SATURATE, else 60000
P.DetThresh     = 8;       % detection threshold [sigma] (all sources, used for isolation)
P.MinSN         = 40;      % minimum S/N for PSF stars
P.MaxPeakFrac   = 0.7;     % reject if peak > MaxPeakFrac*Satur (non-linearity)
P.StampRad      = 10;      % cutout half-size [pix]
P.IsoRad        = 14;      % reject if another source within this radius [pix]
P.IsoFluxFrac   = 0.02;    % ...and that neighbour is brighter than this fraction
P.WinFactor     = 2.5;     % broad moment window sigma = WinFactor * sigma_PSF
P.CoreFactor    = 0.6;     % core window sigma         = CoreFactor * sigma_PSF
P.MapGrid       = [8 12];  % [Nx Ny] cells for maps (portrait, like the LAST frame)
P.PsfGrid       = [4 6];   % [Nx Ny] cells for the ePSF mosaic (portrait, like the LAST frame)
P.Oversamp      = 4;       % ePSF oversampling
P.MaxStarsPsf   = 300;     % max stars per ePSF cell (highest S/N kept)
P.CropGrid      = [8 12];  % [Nx Ny] crops: the image is processed crop by crop to save memory
P.CropMargin    = 32;      % overlap margin [pix] around each crop (>= IsoRad + StampRad)
P.SigmaPSF      = [];      % PSF sigma [pix] for the moment windows; [] -> estimated on the central crop
P.OpticalCenter = [];      % [X Y] pix of the optical axis; [] -> image center
P.UseFlags      = true;    % drop stars with non-zero AstroPack FLAGS
P.FontSize      = 13;      % axes/label font size; titles +3 (bold), legends -1, colorbars -1
P.FigSizeMap    = [760 760];  % [W H] pixels of map figures
P.FigSizePlot   = [820 560];  % [W H] pixels of x-y plot figures
P.MosaicTilePx  = 150;       % tile size [pix] in the ePSF mosaic (figure size follows P.PsfGrid)
P.MosaicTitleFS = 8;         % font size [pt] of the ePSF tile titles
P.BackSubSize   = [128 128];  % background cell size inside each crop

%% ================ Optical design: RASA 11 v2 + QHY600M ================
% Basic first-order parameters (Celestron spec / LAST camera). Verify these
% against your design analysis and replace if different.
D.Aperture_mm    = 279.4;
D.FocalLength_mm = 620;
D.FNumber        = D.FocalLength_mm / D.Aperture_mm;      % ~2.22
D.PixSize_um     = 3.76;                                  % IMX455
D.ImageCircle_mm = 43.3;
D.Lambda_um      = 0.55;
% Predictions from YOUR ray-trace / design analysis (leave [] to skip overlay):
D.Field_mm       = [];   % field radii of the predictions, e.g. [0 5 10 15 20 21.6]
D.RMSspot_um     = [];   % polychromatic RMS spot RADIUS [um] at D.Field_mm
D.ComaShift_um   = [];   % predicted centroid-minus-core shift [um], + = outward

PixScale = 206.265 * D.PixSize_um / D.FocalLength_mm;    % arcsec/pix
CFZ_um   = 4.88 * D.Lambda_um * D.FNumber^2;             % +/- critical focus zone
% --- choose the image file
if isempty(P.ImageFile)
    [FName, FDir] = uigetfile({'*.fits;*.fit;*.fts;*.fits.fz', 'FITS images (*.fits, *.fit, *.fts, *.fits.fz)'; ...
                               '*.*', 'All files'}, 'Select a RASA 11 star-field image', P.ImageDir);
    if isequal(FName, 0); error('No image selected.'); end
    P.File = fullfile(FDir, FName);
else
    if isfile(P.ImageFile)
        P.File = P.ImageFile;                          % full or relative path given
    else
        P.File = fullfile(P.ImageDir, P.ImageFile);    % name inside P.ImageDir
    end
    if ~isfile(P.File); error('Image file not found: %s', P.File); end
end
[ImDir, ImName] = fileparts(P.File);
if isempty(ImDir); ImDir = pwd; end              % bare file name -> current folder
P.OutDir = fullfile(ImDir, [P.OutDirPrefix ImName]);
if ~exist(P.OutDir, 'dir'); mkdir(P.OutDir); end
LogFile = fullfile(P.OutDir, 'console.log');      % console output, embedded in the HTML report
diary off;  if exist(LogFile, 'file'); delete(LogFile); end
diary(LogFile);
RunStart = now;                                  % datenum

%% ============================ Load image ===============================
fprintf('Reading %s\n', P.File);
AI = AstroImage({P.File}, 'HDU', P.HDU);
ImFull = AI.Image;                       % keep native class (e.g. uint16): 4x less memory than double
Satur  = P.Satur;
if isempty(Satur)
    try
        Satur = AI.HeaderData.getVal({'SATURVAL', 'SATURATE'});
    catch
        Satur = NaN;
    end
    if isempty(Satur) || ~isfinite(Satur); Satur = 60000; end
end
clear AI
[Ny, Nx] = size(ImFull);
fprintf('Image %d x %d pix, class %s, saturation %g ADU\n', Nx, Ny, class(ImFull), Satur);

%% ===================== Crop layout & PSF window size ===================
% Interior tiles (non-overlapping) define which stars a crop "owns";
% each crop is read with a margin so stars near tile edges are fully measured
% and the isolation test sees neighbours across the tile boundary.
xe = round(linspace(0, Nx, P.CropGrid(1)+1));
ye = round(linspace(0, Ny, P.CropGrid(2)+1));
NCrop = P.CropGrid(1) * P.CropGrid(2);
R = P.StampRad;

Sig = P.SigmaPSF;
if isempty(Sig)
    % estimate on the central crop: seed from AstroPack 2nd moments, refine twice
    ic = ceil(P.CropGrid(1)/2);  jc = ceil(P.CropGrid(2)/2);
    C  = analyzeCrop(ImFull, xe, ye, ic, jc, P, 'Satur', Satur);
    Sig = C.SigSeed;
    for pass = 1:2
        C = analyzeCrop(ImFull, xe, ye, ic, jc, P, 'Satur', Satur, 'Sig', Sig);
        ok = isfinite(C.S.sigma);
        if nnz(ok) >= 5; Sig = median(C.S.sigma(ok)); end
    end
end
if ~isfinite(Sig) || Sig <= 0; Sig = 1.2; end
fprintf('PSF sigma for moment windows = %.2f pix (FWHM %.2f pix = %.2f")\n', ...
        Sig, 2.3548*Sig, 2.3548*Sig*PixScale);

%% ======================= Crop-by-crop measurement ======================
% ePSF accumulators, one per P.PsfGrid cell
Hs = 6;  nGp = 2*Hs*P.Oversamp + 1;
NPsf = P.PsfGrid(1) * P.PsfGrid(2);
EAcc = zeros(nGp, nGp, NPsf);  ECnt = zeros(nGp, nGp, NPsf);  ENst = zeros(NPsf, 1);
pex = linspace(0, Nx, P.PsfGrid(1)+1);  pey = linspace(0, Ny, P.PsfGrid(2)+1);

Sc = cell(NCrop, 1);
CropRows = zeros(NCrop, 10);   % [k i j x1 x2 y1 y2 NDet NSel t]
NDet = 0;  NSel = 0;  k = 0;  tic
for j = 1:P.CropGrid(2)
    for i = 1:P.CropGrid(1)
        k = k + 1;
        C = analyzeCrop(ImFull, xe, ye, i, j, P, 'Satur', Satur, 'Sig', Sig, ...
                        'Psf', struct('pex', pex, 'pey', pey, 'Hs', Hs, 'Ov', P.Oversamp));
        NDet = NDet + C.NDet;  NSel = NSel + C.NSel;
        Sc{k} = C.S;
        EAcc = EAcc + C.EAcc;  ECnt = ECnt + C.ECnt;  ENst = ENst + C.ENst;
        CropRows(k,:) = [k i j xe(i)+1 xe(i+1) ye(j)+1 ye(j+1) C.NDet C.NSel toc];
        fprintf('crop %3d/%d  [x %5d-%5d, y %5d-%5d]  detected %5d  PSF stars %4d  (%.0f s)\n', ...
                k, NCrop, xe(i)+1, xe(i+1), ye(j)+1, ye(j+1), C.NDet, numel(C.S.sigma), toc);
    end
end
clear C
S  = catStruct(Sc);  clear Sc
fprintf('%d sources detected, %d isolated unsaturated candidates\n', NDet, NSel);

% global quality cuts on the measured quantities
ok = isfinite(S.sigma) & isfinite(S.xcore) & isfinite(S.xbroad) & ...
     S.sigma > 0.5*Sig & S.sigma < 3*Sig & S.flux > 0;
fn = fieldnames(S);
for f = 1:numel(fn); S.(fn{f}) = S.(fn{f})(ok); end
Nst = numel(S.sigma);
fprintf('%d PSF stars after moment quality cuts (median sigma %.2f pix)\n', Nst, median(S.sigma));
if Nst < 50
    warning('Only %d PSF stars: lower P.MinSN or P.DetThresh.', Nst);
end

%% =========================== Field geometry ============================
if isempty(P.OpticalCenter)
    OC = [(Nx+1)/2, (Ny+1)/2];
else
    OC = P.OpticalCenter;
end
mm  = D.PixSize_um / 1000;                       % mm per pixel
xf  = (S.xcore - OC(1)) * mm;                    % field x [mm]
yf  = (S.ycore - OC(2)) * mm;
rf  = hypot(xf, yf);
ux  = xf ./ max(rf, eps);  uy = yf ./ max(rf, eps);
Rmax = hypot(Nx, Ny)/2 * mm;                     % half-diagonal [mm]

% coma centroid shift [um] and its radial/tangential parts
Cx = (S.xbroad - S.xcore) * D.PixSize_um;
Cy = (S.ybroad - S.ycore) * D.PixSize_um;
CRad = Cx.*ux + Cy.*uy;
CTan = -Cx.*uy + Cy.*ux;
cRad = S.c1.*ux + S.c2.*uy;                      % normalized 3rd moment, radial
cTan = -S.c1.*uy + S.c2.*ux;
% spin-2: radial (+) vs tangential (-) elongation
eRad = S.e1.*(ux.^2-uy.^2) + S.e2.*(2*ux.*uy);
FWHM_um = 2.3548 * S.sigma * D.PixSize_um;

%% ============================== Field fits =============================
% --- coma vector field: [Cx;Cy] = C0 + A [x;y]
n  = Nst; z = zeros(n,1); o = ones(n,1);
Ac = [o z xf yf z z; z o z z xf yf];
[pc, wc] = robustLsq(Ac, [Cx; Cy]);
C0   = pc(1:2);                                   % um
Amat = [pc(3) pc(4); pc(5) pc(6)];                % um/mm
kRad  = (Amat(1,1) + Amat(2,2))/2;                % radial coma gradient
kCurl = (Amat(2,1) - Amat(1,2))/2;                % rotational term (should be ~0)
if abs(det(Amat)) > 1e-12
    ComaFree_mm  = -Amat \ C0;
else
    ComaFree_mm  = [NaN; NaN];
end
ComaFree_pix = OC(:) + ComaFree_mm / mm;

% --- same on normalized third moments (independent of centroid method)
[pcn, ~] = robustLsq(Ac, [S.c1; S.c2]);

% --- astigmatism (spin-2)
qn = (xf.^2 - yf.^2)/Rmax^2;  pn = 2*xf.*yf/Rmax^2;
Aa = [o z xf yf z z qn; z o z z xf yf pn];
[pa, ~] = robustLsq(Aa, [S.e1; S.e2]);
E0 = pa(1:2);  kA = pa(7);

% --- focus surface from FWHM (on grid medians to suppress noise)
G = gridStats(S.xcore, S.ycore, [Nx Ny], P.MapGrid, ...
              struct('FWHM', FWHM_um, 'Cx', Cx, 'Cy', Cy, 'e1', S.e1, 'e2', S.e2, ...
                     'c1', S.c1, 'c2', S.c2, 'tr1', S.t1, 'tr2', S.t2));
gxm = (G.xc - OC(1))*mm;  gym = (G.yc - OC(2))*mm;
use = G.n >= 5 & isfinite(G.FWHM);
Ffit = fitFocusSurface(gxm(use), gym(use), G.FWHM(use), G.n(use), D.FNumber);
Tilt_mrad   = hypot(Ffit.tx, Ffit.ty);                     % um/mm == mrad
Tilt_arcmin = Tilt_mrad * 1e-3 * 180/pi * 60;
TiltAxisPA  = atan2d(Ffit.ty, Ffit.tx) + 90;               % line of equal focus
dZ_chip     = Tilt_mrad * hypot(Nx, Ny) * mm;              % um across diagonal
if abs(Ffit.cz) > 1e-9
    Rcurv_mm = 1 / (2*Ffit.cz*1e-3);                        % sag = r^2/(2Rc)
else
    Rcurv_mm = Inf;
end

%% ============================== Summary ================================
% Built once as rows {group, quantity, value, note}; printed to the console
% and reused in the HTML report.
asec = @(um) um/D.PixSize_um*PixScale;            % um on the detector -> arcsec
Info = {'Image file',        P.File;
        'Image size',        sprintf('%d x %d pix', Nx, Ny);
        'Pixel scale',       sprintf('%.3f "/pix  (%.2f um pixels)', PixScale, D.PixSize_um);
        'Focal ratio / CFZ', sprintf('f/%.2f,  CFZ = +/-%.1f um', D.FNumber, CFZ_um);
        'Optical center',    sprintf('(%.1f, %.1f) pix (assumed)', OC);
        'Crop grid',         sprintf('%d x %d crops, %d pix margin', P.CropGrid, P.CropMargin);
        'Sources detected',  sprintf('%d  (%d isolated, unsaturated candidates)', NDet, NSel);
        'PSF stars used',    sprintf('%d', Nst);
        'Moment window PSF', sprintf('sigma %.2f pix, FWHM %.2f"', Sig, 2.3548*Sig*PixScale)};

SumRows = {
 'Coma', 'Constant term C0',     sprintf('(%+.2f, %+.2f) um,  |C0| = %.2f um = %.2f"', C0, norm(C0), asec(norm(C0))), 'collimation / decenter';
 'Coma', 'Radial gradient kRad', sprintf('%+.3f um/mm  ->  %+.2f um at r = %.1f mm', kRad, kRad*Rmax, Rmax), '+ : flare points outward';
 'Coma', 'Rotational term kCurl', sprintf('%+.3f um/mm', kCurl), 'should be ~0';
 'Coma', 'Coma-free point',      sprintf('(%+.2f, %+.2f) mm  =  pix (%.0f, %.0f)', ComaFree_mm, ComaFree_pix), 'where the optical axis lands';
 'Coma', 'Normalized 3rd-moment fit', sprintf('C0 = (%+.3f, %+.3f),  kRad = %+.4f /mm', pcn(1:2), (pcn(3)+pcn(6))/2), 'independent cross-check';
 'Coma', 'Median |coma| corners / center', sprintf('%.2f / %.2f um', ...
        median(hypot(Cx(rf > 0.8*Rmax), Cy(rf > 0.8*Rmax))), median(hypot(Cx(rf < 0.25*Rmax), Cy(rf < 0.25*Rmax)))), 'r > 0.8 Rmax / r < 0.25 Rmax';
 'Astigmatism', 'Constant e0',   sprintf('(%+.3f, %+.3f),  |e0| = %.3f,  PA = %.1f deg', E0, norm(E0), 0.5*atan2d(E0(2), E0(1))), 'tracking / wind / pinching';
 'Astigmatism', 'Field term kA', sprintf('%+.3f at Rmax', kA), '+ radial / - tangential elongation';
 'Focus', 'F0 (seeing + design core)', sprintf('%.2f um = %.2f"', Ffit.F0, asec(Ffit.F0)), '';
 'Focus', 'Defocus at center |z0|',    sprintf('%.1f um', abs(Ffit.z0)), sprintf('CFZ +/-%.1f um', CFZ_um);
 'Focus', 'Sensor tilt',               sprintf('%.2f arcmin,  equal-focus axis PA %.0f deg', Tilt_arcmin, mod(TiltAxisPA,180)), '';
 'Focus', 'Focus change over diagonal', sprintf('%.1f um  (%.1f x CFZ)', dZ_chip, dZ_chip/(2*CFZ_um)), '';
 'Focus', 'Field curvature radius',    sprintf('%.0f mm', Rcurv_mm), '|R| -> inf : flat';
 };

Hints = {};
edgeComa = abs(kRad)*Rmax;
if norm(C0) > 0.5*max(edgeComa, D.PixSize_um/2)
    Hints{end+1} = sprintf(['Constant coma dominates -> collimation/decenter: the optical axis is off the ' ...
        'sensor center by ~%.1f mm (coma-free point). Check camera/corrector tilt and primary ' ...
        'collimation; coma flares point away from the coma-free point.'], norm(ComaFree_mm));
end
if edgeComa > D.PixSize_um/2
    if kRad > 0
        Hints{end+1} = sprintf('Field coma grows outward (flare points outward) by %.1f um at the edge.', edgeComa);
    else
        Hints{end+1} = sprintf('Field coma grows with radius, flare points INWARD (%.1f um at the edge).', edgeComa);
    end
    Hints{end+1} = 'Compare with design; excess field-linear coma suggests a back-focus/spacing error.';
end
if dZ_chip > 2*CFZ_um
    Hints{end+1} = 'Sensor tilt moves focus by more than the CFZ across the chip -> shim the camera.';
end
if abs(kA) > 0.05
    Hints{end+1} = sprintf('Field-dependent astigmatism of %.2f ellipticity at the edge.', kA);
end
if norm(E0) > 0.05
    Hints{end+1} = sprintf('Constant elongation |e0| = %.2f (tracking, wind, or pinched optics).', norm(E0));
end
if isempty(Hints); Hints = {'No aberration exceeds the warning thresholds.'}; end

% per-crop statistics (after the global quality cuts)
CropT = cropStats(CropRows, S.crop, FWHM_um/D.PixSize_um*PixScale, hypot(Cx, Cy), CRad, hypot(S.e1, S.e2));

% console print
fprintf('\n================ RASA 11 aberration summary ================\n');
for r = 1:size(Info,1); fprintf('%-20s: %s\n', Info{r,1}, Info{r,2}); end
grp = '';
for r = 1:size(SumRows,1)
    if ~strcmp(grp, SumRows{r,1}); grp = SumRows{r,1}; fprintf('\n%s\n', upper(grp)); end
    if isempty(SumRows{r,4}); nt = ''; else; nt = ['   [' SumRows{r,4} ']']; end
    fprintf('  %-30s: %s%s\n', SumRows{r,2}, SumRows{r,3}, nt);
end
fprintf('  (the sign of dz is not recoverable from in-focus images)\n');
fprintf('\nDIAGNOSIS HINTS\n');
for h = 1:numel(Hints); fprintf('  * %s\n', Hints{h}); end
fprintf('============================================================\n\n');

%% =============================== Figures ===============================
% One plot per figure. Each is saved as PNG (for the report) and MATLAB .fig.
% FigList collects {base name, caption} in report order.
FigList = cell(0, 2);
rb  = linspace(0, max(rf), 12);

% --- 1. Coma vector map (grid medians) + model + coma-free point
f = newFig('Coma map', P.FigSizeMap, 1);
gc = hypot(G.Cx, G.Cy);
mapImage(gc, Nx, Ny, P.MapGrid);
cb = colorbar; cb.Label.String = '|coma shift| [\mum]';
sc = 0.25 * min(Nx/P.MapGrid(1), Ny/P.MapGrid(2)) / max(gc(isfinite(gc)) + eps);
quiver(G.xc, G.yc, G.Cx*sc, G.Cy*sc, 0, 'k', 'LineWidth', 1.6, 'MaxHeadSize', 0.6);
mx = C0(1) + Amat(1,:)*[gxm gym]';  my = C0(2) + Amat(2,:)*[gxm gym]';
quiver(G.xc, G.yc, mx'*sc, my'*sc, 0, 'r--', 'LineWidth', 1);
plot(OC(1), OC(2), 'w+', 'MarkerSize', 14, 'LineWidth', 2);
plot(ComaFree_pix(1), ComaFree_pix(2), 'mo', 'MarkerSize', 12, 'LineWidth', 2);
lg = legend({'measured (core\rightarrowflare)', 'linear model', 'image center', 'coma-free point'}, ...
       'Location', 'southoutside', 'NumColumns', 2);
xlabel('X [pix]'); ylabel('Y [pix]');
title({'Coma map', sprintf('|C_0| = %.2f \\mum,  k_{rad} = %+.3f \\mum/mm', norm(C0), kRad)}, ...
      'FontSize', P.FontSize, 'FontWeight', 'bold', 'Tag', 'KeepFont', 'UserData', P.FontSize);
FigList(end+1,:) = saveFig(f, P, 'coma_map', ...
    'Coma vector map: grid medians (black, core -> flare), linear field model (red dashed), image center (+) and coma-free point (o).');

% --- 2. Centroid-shift coma vs field radius
f = newFig('Coma vs radius', P.FigSizePlot, 2);  hold on; box on
plot(rf, CRad, '.', 'Color', [0.7 0.7 0.9]);
[rbm, mR] = binMedian(rf, CRad, rb);  [~, mT] = binMedian(rf, CTan, rb);
plot(rbm, mR, 'bo-', 'LineWidth', 1.5); plot(rbm, mT, 'rs-', 'LineWidth', 1.5);
lgs = {'stars (radial)', 'radial median', 'tangential median'};
if ~isempty(D.ComaShift_um)
    plot(D.Field_mm, D.ComaShift_um, 'k--', 'LineWidth', 1.5);  lgs{end+1} = 'design';
end
yline(0, 'k:', 'HandleVisibility', 'off');
legend(lgs, 'Location', 'best');
xlabel('field radius [mm]'); ylabel('coma shift [\mum]');
ylim(qtile([CRad; CTan], [0.02 0.98]));
title('Centroid-shift coma vs field radius');
FigList(end+1,:) = saveFig(f, P, 'coma_vs_radius', ...
    'Centroid-shift coma vs field radius: per-star radial component (dots) and binned medians of the radial and tangential components. Radial > 0 means the flare points outward.');

% --- 3. Normalized third-moment coma vs field radius
f = newFig('Coma moment vs radius', P.FigSizePlot, 3);  hold on; box on
[~, nR] = binMedian(rf, cRad, rb);  [~, nT] = binMedian(rf, cTan, rb);
plot(rbm, nR, 'bo-', 'LineWidth', 1.5); plot(rbm, nT, 'rs-', 'LineWidth', 1.5);
yline(0, 'k:', 'HandleVisibility', 'off');
legend({'radial', 'tangential'}, 'Location', 'best');
xlabel('field radius [mm]'); ylabel('(M_{30}+M_{12}) / s^3');
title('Normalized 3rd-moment coma vs field radius');
FigList(end+1,:) = saveFig(f, P, 'coma_moment_vs_radius', ...
    'Normalized third-moment coma (binned medians), radial and tangential components - an independent cross-check of the centroid-shift measure.');

% --- 4. FWHM map
f = newFig('FWHM map', P.FigSizeMap, 4);
mapImage(G.FWHM/D.PixSize_um*PixScale, Nx, Ny, P.MapGrid);
cb = colorbar; cb.Label.String = 'FWHM ["]';
xlabel('X [pix]'); ylabel('Y [pix]');
title('FWHM map');
FigList(end+1,:) = saveFig(f, P, 'fwhm_map', 'Median FWHM in each map cell.');

% --- 5. FWHM vs field radius
f = newFig('FWHM vs radius', P.FigSizePlot, 5);  hold on; box on
plot(rf, FWHM_um/D.PixSize_um*PixScale, '.', 'Color', [0.75 0.75 0.75]);
[~, mF] = binMedian(rf, FWHM_um, rb);
plot(rbm, mF/D.PixSize_um*PixScale, 'ko-', 'LineWidth', 1.5);
lgs = {'stars', 'binned median'};
if ~isempty(D.RMSspot_um)
    % expected FWHM: seeing core (F0) convolved with design spot (Gaussian approx.)
    Fdes = sqrt(Ffit.F0^2 + (2.3548*D.RMSspot_um/sqrt(2)).^2);
    plot(D.Field_mm, Fdes/D.PixSize_um*PixScale, 'r--', 'LineWidth', 1.5);  lgs{end+1} = 'design + F0';
end
legend(lgs, 'Location', 'best');
xlabel('field radius [mm]'); ylabel('FWHM ["]');
ylim(qtile(FWHM_um/D.PixSize_um*PixScale, [0.01 0.99]));
title('FWHM vs field radius');
FigList(end+1,:) = saveFig(f, P, 'fwhm_vs_radius', 'Per-star FWHM vs field radius with binned medians.');

% --- 6. Focus surface
f = newFig('Focus surface', P.FigSizeMap, 6);
[XXg, YYg] = meshgrid(linspace(-Nx/2, Nx/2, 60)*mm, linspace(-Ny/2, Ny/2, 90)*mm);
Zg = Ffit.z0 + Ffit.tx*XXg + Ffit.ty*YYg + Ffit.cz*(XXg.^2+YYg.^2);
imagesc(XXg(1,:), YYg(:,1), Zg); axis xy image; hold on
contour(XXg, YYg, Zg, [-CFZ_um CFZ_um], 'w--', 'LineWidth', 1.2);
cb = colorbar; cb.Label.String = 'defocus [\mum] (\pm sign ambiguous)';
xlabel('x [mm]'); ylabel('y [mm]');
title(sprintf('Focus surface: tilt %.1f'', \\Deltaz_{diag} = %.0f \\mum', Tilt_arcmin, dZ_chip));
FigList(end+1,:) = saveFig(f, P, 'focus_surface', ...
    'Fitted focus surface (offset + tilt + curvature); white dashed contours mark +/- CFZ.');

% --- 7. Ellipticity whiskers
f = newFig('Ellipticity map', P.FigSizeMap, 7);
ge = hypot(G.e1, G.e2);  gpa = 0.5*atan2(G.e2, G.e1);
mapImage(ge, Nx, Ny, P.MapGrid);
cb = colorbar; cb.Label.String = '|e|';
L  = 0.4 * min(Nx/P.MapGrid(1), Ny/P.MapGrid(2)) / max(ge(isfinite(ge)) + eps);
wx = L*ge.*cos(gpa); wy = L*ge.*sin(gpa);
plot([G.xc-wx G.xc+wx]', [G.yc-wy G.yc+wy]', 'k-', 'LineWidth', 2.5);
xlabel('X [pix]'); ylabel('Y [pix]');
title(sprintf('Ellipticity: |e_0| = %.3f, k_A = %+.3f (+radial)', norm(E0), kA));
FigList(end+1,:) = saveFig(f, P, 'ellipticity_map', 'Ellipticity whiskers (length ~ |e|, direction = major axis).');

% --- 8. ePSF mosaic (accumulated crop by crop) - one figure, one tile per field cell
MosWH = [P.PsfGrid(1)*P.MosaicTilePx + 40, P.PsfGrid(2)*(P.MosaicTilePx + 14) + 60];
f = newFig('ePSF mosaic', MosWH, 8);
try
    tiledlayout(P.PsfGrid(2), P.PsfGrid(1), 'TileSpacing', 'tight', 'Padding', 'tight');      % R2021a+
catch
    tiledlayout(P.PsfGrid(2), P.PsfGrid(1), 'TileSpacing', 'compact', 'Padding', 'compact');
end
ax = (-Hs*P.Oversamp:Hs*P.Oversamp)/P.Oversamp;
for iy = P.PsfGrid(2):-1:1                   % top row first = high Y
    for ix = 1:P.PsfGrid(1)
        nexttile;
        c  = (iy-1)*P.PsfGrid(1) + ix;
        E  = EAcc(:,:,c) ./ max(ECnt(:,:,c), 1);
        E  = conv2(E, ones(P.Oversamp)/P.Oversamp^2, 'same');   % fill fine-grid holes
        in = S.xcore > pex(ix) & S.xcore <= pex(ix+1) & S.ycore > pey(iy) & S.ycore <= pey(iy+1);
        imagesc(ax, ax, log10(max(E/max(E(:)), 1e-3))); axis xy image; colormap(gca, 'parula');
        hold on; contour(ax, ax, E/max(E(:)), [0.5 0.2 0.05 0.02], 'w', 'LineWidth', 0.6);
        mcx = median(Cx(in))/D.PixSize_um;  mcy = median(Cy(in))/D.PixSize_um;
        quiver(0, 0, mcx*8, mcy*8, 0, 'r', 'LineWidth', 1.8, 'MaxHeadSize', 1);
        title(sprintf('n=%d, FWHM %.1f"', ENst(c), median(FWHM_um(in))/D.PixSize_um*PixScale), ...
              'FontSize', P.MosaicTitleFS, 'FontWeight', 'normal', 'Tag', 'KeepFont', 'UserData', P.MosaicTitleFS);
        set(gca, 'XTick', [], 'YTick', []);
    end
end
sgtitle('Oversampled ePSF across the field', 'FontSize', P.FontSize+1, 'FontWeight', 'bold');
FigList(end+1,:) = saveFig(f, P, 'epsf_mosaic', ...
    'Oversampled ePSF in each field cell (log scale), contours at 50/20/5/2% of peak; red arrow = median coma x8.');

%% =============================== Save ==================================
T = table(S.crop, S.xcore, S.ycore, xf, yf, rf, S.sn, S.flux, ...
          FWHM_um/D.PixSize_um*PixScale, S.e1, S.e2, eRad, Cx, Cy, CRad, CTan, ...
          S.c1, S.c2, cRad, cTan, S.t1, S.t2, ...
          'VariableNames', {'Crop','X','Y','x_mm','y_mm','r_mm','SN','Flux', ...
          'FWHM_arcsec','e1','e2','eRad','ComaX_um','ComaY_um','ComaRad_um','ComaTan_um', ...
          'c1','c2','cRad','cTan','Trefoil1','Trefoil2'});
writetable(T, fullfile(P.OutDir, 'psf_stars.csv'));
Res = struct('P', P, 'D', D, 'PixScale', PixScale, 'OpticalCenter', OC, ...
             'C0_um', C0, 'A_um_per_mm', Amat, 'kRad', kRad, 'kCurl', kCurl, ...
             'ComaFree_mm', ComaFree_mm, 'ComaFree_pix', ComaFree_pix, ...
             'NormComaFit', pcn, 'E0', E0, 'kA', kA, 'AstigFit', pa, ...
             'Focus', Ffit, 'Tilt_arcmin', Tilt_arcmin, 'dZ_diag_um', dZ_chip, ...
             'Rcurv_mm', Rcurv_mm, 'Grid', G, 'SigmaPSF', Sig, ...
             'ePSF_Acc', EAcc, 'ePSF_Cnt', ECnt, 'ePSF_N', ENst);
[~, ImBase, ImExt] = fileparts(P.File);
MatFile = fullfile(P.OutDir, [lastBaseName([ImBase ImExt]) '_IQA.mat']);
save(MatFile, 'Res', 'T');
fprintf('Results written to %s\n', P.OutDir);
diary off;

%% ============================ HTML report ==============================
Figs = FigList;
RepFile = fullfile(P.OutDir, [lastBaseName([ImBase ImExt]) '_report.html']);
writeReport(RepFile, P, Info, SumRows, Hints, CropT, Figs, LogFile, RunStart);
fprintf('HTML report: %s\n', RepFile);


%% ============================ Local functions ==========================
function T = cropStats(CropRows, crop, fwhm, coma, comaRad, ell)
% Per-crop summary (struct of column vectors): layout, detection counts, PSF-star medians.
    n = size(CropRows, 1);
    med = @(v) accumarray(crop, v, [n 1], @(x) median(x, 'omitnan'), NaN);
    T = struct('Crop', CropRows(:,1), 'Col', CropRows(:,2), 'Row', CropRows(:,3), ...
               'X1', CropRows(:,4), 'X2', CropRows(:,5), 'Y1', CropRows(:,6), 'Y2', CropRows(:,7), ...
               'Detected', CropRows(:,8), 'Candidates', CropRows(:,9), ...
               'PSFstars', accumarray(crop, 1, [n 1]), ...
               'FWHM_arcsec', med(fwhm), 'Coma_um', med(coma), 'ComaRad_um', med(comaRad), ...
               'Ellip', med(ell), 't_s', CropRows(:,10));
end

function writeReport(File, P, Info, SumRows, Hints, CropT, Figs, LogFile, RunStart)
% Self-contained HTML report (figures embedded as base64 PNG).
    fid = fopen(File, 'w', 'n', 'UTF-8');
    if fid < 0; warning('Cannot write %s', File); return; end
    c = onCleanup(@() fclose(fid));
    w = @(varargin) fprintf(fid, varargin{:});
    [~, nm, ext] = fileparts(P.File);
    Title = reportTitle([nm ext]);

    w('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">\n');
    w('<meta name="viewport" content="width=device-width, initial-scale=1">\n');
    w('<title>%s</title>\n', esc(Title));
    w('%s', reportCss());
    w('</head><body><main>\n');
    w('<header><h1>%s</h1>\n', esc(Title));
    w('<p class="sub">%s &middot; generated %s</p></header>\n', esc([nm ext]), datestr(now, 'yyyy-mm-dd HH:MM'));

    % --- run info
    w('<section><h2>Image &amp; run</h2><table class="kv">\n');
    for r = 1:size(Info,1)
        w('<tr><th>%s</th><td>%s</td></tr>\n', esc(Info{r,1}), units(esc(Info{r,2})));
    end
    w('<tr><th>Run time</th><td>%.0f s</td></tr>\n', (now - RunStart)*86400);
    w('</table></section>\n');

    % --- diagnosis
    w('<section><h2>Diagnosis</h2><ul class="hints">\n');
    for h = 1:numel(Hints); w('<li>%s</li>\n', units(esc(Hints{h}))); end
    w('</ul></section>\n');

    % --- aberration summary
    w('<section><h2>Aberration summary</h2><table class="sum"><thead><tr><th>Quantity</th><th>Value</th><th>Note</th></tr></thead><tbody>\n');
    grp = '';
    for r = 1:size(SumRows,1)
        if ~strcmp(grp, SumRows{r,1})
            grp = SumRows{r,1};
            w('<tr class="grp"><td colspan="3">%s</td></tr>\n', esc(grp));
        end
        w('<tr><td>%s</td><td class="num">%s</td><td class="note">%s</td></tr>\n', ...
          esc(SumRows{r,2}), units(esc(SumRows{r,3})), units(esc(SumRows{r,4})));
    end
    w('</tbody></table><p class="foot">The sign of defocus is not recoverable from in-focus images.</p></section>\n');

    % --- figures: one per row
    w('<section><h2>Figures</h2><p class="foot">Each figure is also saved as a MATLAB .fig file in the report folder.</p>\n<div class="figs">\n');
    for f = 1:size(Figs,1)
        fn = fullfile(P.OutDir, [Figs{f,1} '.png']);
        if ~exist(fn, 'file'); continue; end
        fi = fopen(fn, 'r');  b = fread(fi, Inf, 'uint8=>uint8');  fclose(fi);
        w('<figure><img alt="%s" src="data:image/png;base64,%s"><figcaption><b>Fig. %d.</b> %s ', ...
          esc(Figs{f,1}), matlab.net.base64encode(b), f, units(esc(Figs{f,2})));
        w('<a href="%s.fig">%s.fig</a></figcaption></figure>\n', Figs{f,1}, Figs{f,1});
    end
    w('</div></section>\n');

    % --- crop table
    w('<section><h2>Crops</h2>\n<p class="foot">Click a column header to sort. Shading: darker = larger value in that column. ');
    w('Medians are over the PSF stars of each crop after the global quality cuts.</p>\n');
    w('<div class="scroll"><table id="crops" class="crops"><thead><tr>');
    hdr = {'Crop','Col','Row','X range','Y range','Detected','Candidates','PSF stars', ...
           'FWHM ["]','|Coma| [&micro;m]','Radial coma [&micro;m]','|e|','t [s]'};
    for h = 1:numel(hdr); w('<th data-col="%d">%s</th>', h-1, hdr{h}); end
    w('</tr></thead><tbody>\n');
    shade = {'FWHM_arcsec', 'Coma_um', 'Ellip'};
    rng = struct();
    for q = 1:numel(shade)
        v = CropT.(shade{q});  rng.(shade{q}) = [min(v(isfinite(v))) max(v(isfinite(v)))];
        if isempty(rng.(shade{q})); rng.(shade{q}) = [0 1]; end
    end
    fnC = fieldnames(CropT);
    for r = 1:numel(CropT.Crop)
        for q = 1:numel(fnC); t.(fnC{q}) = CropT.(fnC{q})(r); end
        w('<tr><td>%d</td><td>%d</td><td>%d</td><td data-v="%d">%d&ndash;%d</td><td data-v="%d">%d&ndash;%d</td>', ...
          t.Crop, t.Col, t.Row, t.X1, t.X1, t.X2, t.Y1, t.Y1, t.Y2);
        w('<td>%d</td><td>%d</td><td>%d</td>', t.Detected, t.Candidates, t.PSFstars);
        w('%s%s%s%s', heatCell(t.FWHM_arcsec, rng.FWHM_arcsec, '%.2f'), heatCell(t.Coma_um, rng.Coma_um, '%.2f'), ...
          plainCell(t.ComaRad_um, '%+.2f'), heatCell(t.Ellip, rng.Ellip, '%.3f'));
        w('<td>%.0f</td></tr>\n', t.t_s);
    end
    w('</tbody></table></div></section>\n');

    % --- console log
    w('<section><details><summary>Console output</summary><pre>');
    if exist(LogFile, 'file'); w('%s', esc(fileread(LogFile))); end
    w('</pre></details></section>\n');

    w('</main>\n%s\n</body></html>\n', reportJs());
end

function mapImage(V, Nx, Ny, Grid)
% Draw grid-cell values so the cells exactly tile the detector, with the
% axes limited to the image pixels [1..Nx] x [1..Ny].
    dx = Nx/Grid(1);  dy = Ny/Grid(2);
    imagesc([dx/2, Nx-dx/2], [dy/2, Ny-dy/2], reshape(V, Grid(2), Grid(1)));
    axis xy image;  hold on
    xlim([1 Nx]);  ylim([1 Ny]);
    set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');   % keep limits when arrows/markers are added
end

function f = newFig(Name, WH, k)
% New figure with a cascaded position and the requested size [W H] in pixels.
    f = figure('Name', Name, 'Color', 'w', 'Position', [40+25*k, 40+20*k, WH(1), WH(2)]);
end

function Item = saveFig(f, P, Base, Caption)
% Apply uniform fonts, then save PNG + MATLAB .fig. Returns {Base, Caption}.
    styleFig(f, P.FontSize);
    exportgraphics(f, fullfile(P.OutDir, [Base '.png']), 'Resolution', 150);
    savefig(f, fullfile(P.OutDir, [Base '.fig']));
    Item = {Base, Caption};
end

function styleFig(f, FS)
% Uniform fonts: axes/labels FS, titles FS+3 bold, legends and colorbars FS-1.
% Titles tagged 'KeepFont' keep their own size, stored in the title's UserData
% (coma map, ePSF tiles). It must be re-applied AFTER setting the axes
% FontSize, because MATLAB rescales the title to FontSize*TitleFontSizeMultiplier.
% Legends that stick out of the figure are moved back inside.
    for a = findobj(f, 'Type', 'axes')'
        ht  = get(a, 'Title');
        keep = strcmp(get(ht, 'Tag'), 'KeepFont');
        set(a, 'FontSize', FS, 'LineWidth', 0.8);
        set([get(a, 'XLabel') get(a, 'YLabel')], 'FontSize', FS);
        if keep
            tfs = get(ht, 'UserData');
            if isnumeric(tfs) && isscalar(tfs); set(ht, 'FontSize', tfs); end
        else
            set(ht, 'FontSize', FS + 3, 'FontWeight', 'bold');
        end
    end
    set(findobj(f, 'Type', 'legend'),   'FontSize', FS - 1);
    set(findobj(f, 'Type', 'colorbar'), 'FontSize', FS - 1);
    drawnow;
    for lg = findobj(f, 'Type', 'legend')'
        set(lg, 'Units', 'normalized');
        lp = get(lg, 'Position');
        nx = min(max(lp(1), 0.01), 0.99 - lp(3));          % clamp into [0.01, 0.99]
        if abs(nx - lp(1)) > 1e-3
            for ax = findobj(f, 'Type', 'axes')'           % freeze axes before detaching legend
                set(ax, 'Position', get(ax, 'Position'));
            end
            set(lg, 'Position', [nx lp(2:4)]);
        end
    end
end

function Base = lastBaseName(FileName)
% Base name for output files, built from the same fields as reportTitle:
% 'LAST.01.<mount>.<tel>_yyyymmdd...fits' -> 'LAST_<mount>.<tel>_yyyy-mm-dd'.
% Other names -> the image file name without extension.
% Used for '<Base>_report.html' and '<Base>_IQA.mat'.
    tok = regexp(FileName, '^LAST\.\d+\.(\d+)\.(\d+)_(\d{4})(\d{2})(\d{2})', 'tokens', 'once');
    if isempty(tok)
        [~, Base] = fileparts(FileName);
    else
        Base = sprintf('LAST_%s.%s_%s-%s-%s', tok{:});
    end
end

function T = reportTitle(FileName)
% 'LAST.01.<mount>.<tel>_yyyymmdd...fits' -> 'LAST <mount>.<tel> Telescope Report for yyyy-mm-dd'
    tok = regexp(FileName, '^LAST\.\d+\.(\d+)\.(\d+)_(\d{4})(\d{2})(\d{2})', 'tokens', 'once');
    if isempty(tok)
        T = sprintf('RASA 11 Aberration Report - %s', FileName);   % name not in LAST format
    else
        T = sprintf('LAST %s.%s Telescope Report for %s-%s-%s', tok{:});
    end
end

function s = heatCell(v, rng, fmt)
    if ~isfinite(v); s = '<td class="na">&ndash;</td>'; return; end
    a = 0.06 + 0.5 * (v - rng(1)) / max(rng(2) - rng(1), eps);
    s = sprintf(['<td data-v="%.6g" style="background:rgba(217,95,2,%.2f)">' fmt '</td>'], v, a, v);
end

function s = plainCell(v, fmt)
    if ~isfinite(v); s = '<td class="na">&ndash;</td>'; return; end
    s = sprintf(['<td data-v="%.6g">' fmt '</td>'], v, v);
end

function s = esc(s)
    s = strrep(strrep(strrep(char(s), '&', '&amp;'), '<', '&lt;'), '>', '&gt;');
end

function s = units(s)
% "um" -> micron sign, "+/-" -> plus-minus, "->" -> arrow (on escaped text)
    s = regexprep(s, '\<um\>', '&micro;m');
    s = strrep(s, '+/-', '&plusmn;');
    s = strrep(s, '-&gt;', '&rarr;');
end

function s = reportCss()
    s = ['<style>', ...
    ':root{--bg:#fbfbfa;--fg:#1d1f21;--mut:#6b6f75;--line:#e2e2df;--card:#fff;--acc:#d95f02;--grp:#f1f0ec}', ...
    '@media (prefers-color-scheme:dark){:root{--bg:#16181b;--fg:#e6e6e3;--mut:#9aa0a6;--line:#2c3036;--card:#1d2024;--acc:#f28e2b;--grp:#24282d}}', ...
    '*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);', ...
    'font:15px/1.5 -apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}', ...
    'main{max-width:1100px;margin:0 auto;padding:24px 16px 64px}', ...
    'h1{font-size:26px;margin:0}h2{font-size:18px;margin:0 0 12px;padding-bottom:6px;border-bottom:2px solid var(--acc);display:inline-block}', ...
    '.sub,.foot{color:var(--mut);font-size:13px}header{margin-bottom:24px}', ...
    'section{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:18px 20px;margin:0 0 20px}', ...
    'table{border-collapse:collapse;width:100%}th,td{padding:6px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}', ...
    '.kv th{width:200px;color:var(--mut);font-weight:500}', ...
    '.sum .grp td{background:var(--grp);font-weight:600;letter-spacing:.02em}.num{font-variant-numeric:tabular-nums;white-space:nowrap}', ...
    '.note{color:var(--mut);font-size:13px}', ...
    '.hints{margin:0;padding-left:20px}.hints li{margin:4px 0}', ...
    '.figs{display:flex;flex-direction:column;gap:28px}', ...
    'figure{margin:0}figure img{display:block;margin:0 auto;max-width:min(100%,860px);height:auto;border:1px solid var(--line);border-radius:6px;background:#fff}', ...
    'figcaption a{color:var(--acc)}', ...
    'figcaption{color:var(--mut);font-size:inherit;margin-top:8px}', ...
    '.scroll{overflow-x:auto}.crops{font-size:13px;font-variant-numeric:tabular-nums}', ...
    '.crops th{cursor:pointer;white-space:nowrap;position:sticky;top:0;background:var(--card)}.crops th:hover{color:var(--acc)}', ...
    '.crops td{white-space:nowrap}.na{color:var(--mut)}', ...
    'pre{font:12px/1.4 ui-monospace,Menlo,Consolas,monospace;white-space:pre-wrap;max-height:480px;overflow:auto}', ...
    'summary{cursor:pointer;font-weight:600}', ...
    '</style>'];
end

function s = reportJs()
    s = ['<script>(function(){var t=document.getElementById("crops");if(!t)return;', ...
    'var dir={};t.querySelectorAll("th").forEach(function(th){th.addEventListener("click",function(){', ...
    'var c=+th.dataset.col,b=t.tBodies[0],rows=[].slice.call(b.rows);dir[c]=-(dir[c]||-1);', ...
    'var val=function(r){var d=r.cells[c],v=d.dataset.v!==undefined?d.dataset.v:d.textContent;var n=parseFloat(v);return isNaN(n)?-Infinity:n;};', ...
    'rows.sort(function(a,b2){return dir[c]*(val(a)-val(b2));});rows.forEach(function(r){b.appendChild(r);});});});})();</script>'];
end

function S = measureAll(Sub, XY, Args)
    % Window-weighted moments (centroids, size, ellipticity, coma, trefoil) of stars.
    %   For each star a square stamp of half-size StampRad is cut around the
    %   rounded input position and two iterative Gaussian-windowed centroids
    %   are computed (winCentroid): a broad one (window sigma SigWin) and a
    %   core one (window sigma SigCore). Their difference is the coma
    %   centroid shift. Second moments are measured with the broad window
    %   and the window is removed analytically (exact for a Gaussian PSF);
    %   third moments are normalized by the windowed size s^3.
    % Input  : - Background-subtracted image (2D, double). Rows = Y, columns = X.
    %          - [X, Y] matrix (N x 2) of initial star positions [pix], in the
    %            pixel coordinates of the first input. Each position must be
    %            at least StampRad+1 pixels from the image edges.
    %          * ...,key,val,...
    %            'StampRad' - Stamp half-size [pix]; the stamp is (2R+1)^2.
    %                   Default is 10.
    %            'SigWin' - Sigma [pix] of the broad Gaussian window, used for
    %                   the broad centroid and for all 2nd/3rd moments.
    %                   Recommended ~2.5 x PSF sigma. Default is 3.
    %            'SigCore' - Sigma [pix] of the narrow Gaussian window used for
    %                   the core centroid. Recommended ~0.6 x PSF sigma.
    %                   Default is 0.7.
    %            'MaxIter' - Maximum centroid iterations. Default is 30.
    %            'MaxShift' - A star is rejected (NaN) if its broad centroid
    %                   moved more than this from the rounded input position
    %                   [pix]. Default is StampRad/2 (when empty).
    % Output : - Structure of N x 1 column vectors (NaN for rejected stars):
    %            .xcore, .ycore   - core (narrow-window) centroid [pix]
    %            .xbroad, .ybroad - broad-window centroid [pix]
    %            .sigma  - window-deconvolved PSF sigma, det(C)^(1/4) [pix]
    %            .e1, .e2 - ellipticity (Cxx-Cyy)/(Cxx+Cyy), 2Cxy/(Cxx+Cyy)
    %            .c1, .c2 - normalized spin-1 (coma) moments
    %                       (M30+M12)/s^3, (M21+M03)/s^3
    %            .t1, .t2 - normalized spin-3 (trefoil) moments
    %                       (M30-3M12)/s^3, (3M21-M03)/s^3
    %            .flux   - sum of the stamp within radius StampRad [ADU]
    % Example: S = measureAll(Sub, [512.3 400.8; 120.1 77.5], 'StampRad', 10, ...
    %                         'SigWin', 3.5, 'SigCore', 0.85);
    arguments
        Sub
        XY(:,2)
        Args.StampRad(1,1) = 10;
        Args.SigWin(1,1)   = 3;
        Args.SigCore(1,1)  = 0.7;
        Args.MaxIter(1,1)  = 30;
        Args.MaxShift      = [];
    end
    R = Args.StampRad;
    if isempty(Args.MaxShift); MaxShift = R/2; else; MaxShift = Args.MaxShift; end
    [gx, gy] = meshgrid(-R:R, -R:R);
    sigW = Args.SigWin;
    N = size(XY, 1);
    nanv = NaN(N,1);
    S = struct('xcore',nanv,'ycore',nanv,'xbroad',nanv,'ybroad',nanv, ...
               'sigma',nanv,'e1',nanv,'e2',nanv,'c1',nanv,'c2',nanv, ...
               't1',nanv,'t2',nanv,'flux',nanv);
    for i = 1:N
        xi = round(XY(i,1)); yi = round(XY(i,2));
        St = Sub(yi-R:yi+R, xi-R:xi+R);          % rows = Y, cols = X
        X  = gx + xi;  Y = gy + yi;
        [xc, yc] = winCentroid(St, X, Y, XY(i,1), XY(i,2), 'Sigma', sigW,         'MaxIter', Args.MaxIter);
        [xk, yk] = winCentroid(St, X, Y, xc,      yc,      'Sigma', Args.SigCore, 'MaxIter', Args.MaxIter);
        if ~isfinite(xc) || ~isfinite(xk) || hypot(xc-xi, yc-yi) > MaxShift; continue; end
        dx = X - xc;  dy = Y - yc;
        W  = exp(-(dx.^2 + dy.^2)/(2*sigW^2)) .* St;
        F  = sum(W, 'all');  if F <= 0; continue; end
        mxx = sum(W.*dx.^2,'all')/F;  myy = sum(W.*dy.^2,'all')/F;  mxy = sum(W.*dx.*dy,'all')/F;
        m30 = sum(W.*dx.^3,'all')/F;  m03 = sum(W.*dy.^3,'all')/F;
        m21 = sum(W.*dx.^2.*dy,'all')/F;  m12 = sum(W.*dx.*dy.^2,'all')/F;
        Cm = [mxx mxy; mxy myy];
        if det(Cm) <= 0; continue; end
        Ct = inv(inv(Cm) - eye(2)/sigW^2);       % remove Gaussian window (exact for Gaussian PSF)
        if any(eig(Ct) <= 0); continue; end
        s  = sqrt((mxx + myy)/2);                % windowed size for 3rd-moment normalization
        S.xcore(i) = xk;  S.ycore(i) = yk;  S.xbroad(i) = xc;  S.ybroad(i) = yc;
        S.sigma(i) = det(Ct)^(1/4);
        S.e1(i) = (Ct(1,1)-Ct(2,2))/trace(Ct);  S.e2(i) = 2*Ct(1,2)/trace(Ct);
        S.c1(i) = (m30 + m12)/s^3;  S.c2(i) = (m21 + m03)/s^3;
        S.t1(i) = (m30 - 3*m12)/s^3;  S.t2(i) = (3*m21 - m03)/s^3;
        S.flux(i) = sum(St(gx.^2 + gy.^2 <= R^2), 'all');
    end
end

function [Xc, Yc] = winCentroid(Stamp, X, Y, X0, Y0, Args)
    % Iterative Gaussian-windowed centroid of a stamp.
    %   At each iteration the stamp is multiplied by a circular Gaussian of
    %   width Sigma centred on the current estimate, and the intensity-
    %   weighted mean position becomes the new estimate. Iterates until the
    %   step is below Tol or MaxIter is reached.
    % Input  : - Stamp (2D matrix, background subtracted).
    %          - Matrix (same size as the stamp) of the X coordinate [pix] of
    %            each stamp pixel (e.g., from meshgrid).
    %          - Matrix (same size as the stamp) of the Y coordinate [pix] of
    %            each stamp pixel.
    %          - Initial X guess [pix], in the coordinates of the 2nd input.
    %          - Initial Y guess [pix].
    %          * ...,key,val,...
    %            'Sigma' - Sigma of the Gaussian window [pix]. Default is 1.5.
    %            'MaxIter' - Maximum number of iterations. Default is 30.
    %            'MaxStep' - Maximum step per iteration [pix]; larger steps
    %                   are scaled down to this length. Default is 1.
    %            'Tol' - Convergence tolerance on the step length [pix].
    %                   Default is 1e-4.
    % Output : - X centroid [pix] (NaN if the windowed flux is <= 0).
    %          - Y centroid [pix] (NaN if the windowed flux is <= 0).
    % Example: [gx, gy] = meshgrid(-10:10);  St = exp(-((gx-0.3).^2+gy.^2)/2);
    %          [Xc, Yc] = winCentroid(St, gx, gy, 0, 0, 'Sigma', 2);
    arguments
        Stamp
        X
        Y
        X0(1,1)
        Y0(1,1)
        Args.Sigma(1,1)   = 1.5;
        Args.MaxIter(1,1) = 30;
        Args.MaxStep(1,1) = 1;
        Args.Tol(1,1)     = 1e-4;
    end
    Xc = X0;  Yc = Y0;  sig2 = 2*Args.Sigma^2;
    for it = 1:Args.MaxIter
        W = exp(-((X-Xc).^2 + (Y-Yc).^2)/sig2) .* Stamp;
        F = sum(W, 'all');
        if F <= 0; Xc = NaN; Yc = NaN; return; end
        xn = sum(W.*X, 'all')/F;  yn = sum(W.*Y, 'all')/F;
        d  = hypot(xn-Xc, yn-Yc);
        if d > Args.MaxStep                       % limit step
            xn = Xc + (xn-Xc)*Args.MaxStep/d;  yn = Yc + (yn-Yc)*Args.MaxStep/d;
        end
        Xc = xn; Yc = yn;
        if d < Args.Tol; break; end
    end
end

function [p, w] = robustLsq(A, b, nIter, kclip)
% Least squares with iterative 3-sigma (MAD) clipping.
    if nargin < 3; nIter = 6; end
    if nargin < 4; kclip = 3; end
    w = all(isfinite([A b]), 2);
    for it = 1:nIter
        p = A(w,:) \ b(w);
        r = b - A*p;
        s = 1.4826*median(abs(r(w) - median(r(w))));
        wn = isfinite(r) & abs(r) < kclip*s;
        if isequal(wn, w); break; end
        w = wn;
    end
end

function G = gridStats(x, y, NxNy, Grid, V)
% Median of each field in V over a regular grid. Output vectors are ordered
% so that reshape(v, Grid(2), Grid(1)) gives an image (rows = Y).
    ex = linspace(0, NxNy(1), Grid(1)+1);  ey = linspace(0, NxNy(2), Grid(2)+1);
    [IX, IY] = meshgrid(1:Grid(1), 1:Grid(2));
    nc = numel(IX);
    G.xc = (ex(IX(:)) + ex(IX(:)+1))'/2;  G.yc = (ey(IY(:)) + ey(IY(:)+1))'/2;
    G.xc = G.xc(:); G.yc = G.yc(:);  G.n = zeros(nc,1);
    fn = fieldnames(V);
    for f = 1:numel(fn); G.(fn{f}) = NaN(nc,1); end
    for c = 1:nc
        in = x > ex(IX(c)) & x <= ex(IX(c)+1) & y > ey(IY(c)) & y <= ey(IY(c)+1);
        G.n(c) = nnz(in);
        if G.n(c) < 3; continue; end
        for f = 1:numel(fn); G.(fn{f})(c) = median(V.(fn{f})(in), 'omitnan'); end
    end
end

function q = qtile(x, p)
% Quantiles without the Statistics Toolbox.
    x = sort(x(isfinite(x)));  n = numel(x);
    q = interp1(linspace(0, 1, n), x, p);
    if q(2) <= q(1); q(2) = q(1) + 1; end
end

function [xm, ym] = binMedian(x, y, edges)
    xm = NaN(numel(edges)-1,1);  ym = xm;
    for k = 1:numel(edges)-1
        in = x >= edges(k) & x < edges(k+1) & isfinite(y);
        if nnz(in) >= 5; xm(k) = median(x(in)); ym(k) = median(y(in)); end
    end
end

function F = fitFocusSurface(x, y, fwhm, n, N)
% FWHM^2 = F0^2 + ((z0 + tx x + ty y + cz r^2)/N)^2, multi-start fminsearch.
    r2 = x.^2 + y.^2;  w = sqrt(n);
    model = @(p) sqrt(p(1)^2 + ((p(2) + p(3)*x + p(4)*y + p(5)*r2)/N).^2);
    cost  = @(p) sum(w.*(fwhm - model(p)).^2);
    best = Inf;  opt = optimset('MaxFunEvals', 2e4, 'MaxIter', 2e4, 'TolX', 1e-6, 'TolFun', 1e-8);
    F0 = min(fwhm);
    starts = [F0 10 0.5 0.5 0.02; F0 30 -0.5 0.5 -0.02; F0 5 0.5 -0.5 0.05; ...
              0.8*F0 20 1 0 0; 0.8*F0 20 0 1 0; F0 1 0.1 0.1 0.1];
    for s = 1:size(starts,1)
        [p, c] = fminsearch(cost, starts(s,:), opt);
        if c < best; best = c; pb = p; end
    end
    F = struct('F0', abs(pb(1)), 'z0', pb(2), 'tx', pb(3), 'ty', pb(4), 'cz', pb(5), ...
               'rms_um', sqrt(best/sum(w)));
end

function C = analyzeCrop(ImFull, Xe, Ye, I, J, P, Args)
    % Detect, select and measure isolated PSF stars in one crop of a large image.
    %   The crop is the interior tile (Xe(I), Xe(I+1)] x (Ye(J), Ye(J+1)]
    %   extended by P.CropMargin pixels on every side (clipped at the image
    %   edges). The crop gets its own AstroPack background
    %   (imProc.background.background) and source detection
    %   (imProc.sources.findMeasureSources). Candidates are sources whose
    %   position lies inside the interior tile (so each star is owned by
    %   exactly one crop), with S/N >= P.MinSN, no FLAGS, isolated within
    %   P.IsoRad (neighbours searched in the full crop including the margin)
    %   and peak < P.MaxPeakFrac * Satur. Their moments are measured with
    %   measureAll, and optionally added to the ePSF accumulators (accumEPSF).
    % Input  : - Full image (2D, any numeric class, e.g. uint16). Only the
    %            crop is converted to double.
    %          - Vector of tile edges in X [pix], 0..Nx, length Ncol+1.
    %          - Vector of tile edges in Y [pix], 0..Ny, length Nrow+1.
    %          - Tile column index I (1..Ncol).
    %          - Tile row index J (1..Nrow).
    %          - Parameter structure of the main script. Fields used:
    %            StampRad, CropMargin, BackSubSize, DetThresh, MinSN,
    %            UseFlags, IsoRad, IsoFluxFrac, MaxPeakFrac, WinFactor,
    %            CoreFactor, MaxStarsPsf, CropGrid.
    %          * ...,key,val,...
    %            'Satur' - Saturation level [ADU]. Default is 60000.
    %            'Sig' - PSF sigma [pix] used to scale the moment windows
    %                   (SigWin = Sig*P.WinFactor, SigCore = Sig*P.CoreFactor).
    %                   If empty, the moments are not measured and only
    %                   C.SigSeed is returned. Default is [].
    %            'Psf' - ePSF accumulation settings, a structure with fields
    %                   .pex, .pey - ePSF cell edges in X and Y [pix] of the
    %                                full image,
    %                   .Hs        - ePSF half-size [pix],
    %                   .Ov        - oversampling factor.
    %                   If empty, no ePSF is accumulated. Default is [].
    % Output : - Structure with fields:
    %            .NDet    - number of detected sources inside the interior tile.
    %            .NSel    - number of isolated, unsaturated candidates.
    %            .SigSeed - median sigma [pix] from the AstroPack X2/Y2
    %                       moments of the candidates (seed for 'Sig').
    %            .S       - star structure from measureAll, in full-image
    %                       coordinates, plus .sn (S/N, SN_4) and .crop
    %                       (crop index (J-1)*Ncol + I). Empty if 'Sig' is [].
    %            .EAcc    - ePSF sums, nG x nG x Ncell (0 if 'Psf' is []).
    %            .ECnt    - ePSF sample counts, same size as .EAcc.
    %            .ENst    - Ncell x 1 number of stars added per ePSF cell.
    % Example: C = analyzeCrop(ImFull, xe, ye, 4, 6, P, 'Satur', 60000, 'Sig', 1.4, ...
    %                          'Psf', struct('pex', pex, 'pey', pey, 'Hs', 6, 'Ov', 4));
    arguments
        ImFull
        Xe
        Ye
        I(1,1)
        J(1,1)
        P struct
        Args.Satur(1,1) = 60000;
        Args.Sig        = [];
        Args.Psf        = [];
    end
    xe = Xe;  ye = Ye;  i = I;  j = J;
    Satur = Args.Satur;  Sig = Args.Sig;  Psf = Args.Psf;
    [Ny, Nx] = size(ImFull);
    R  = P.StampRad;  M = P.CropMargin;
    x1 = max(1, xe(i)+1-M);  x2 = min(Nx, xe(i+1)+M);
    y1 = max(1, ye(j)+1-M);  y2 = min(Ny, ye(j+1)+M);
    Im = double(ImFull(y1:y2, x1:x2));
    [ny, nx] = size(Im);

    AI = AstroImage({Im});
    AI = imProc.background.background(AI, 'SubSizeXY', min(P.BackSubSize, [nx ny]));
    AI = imProc.sources.findMeasureSources(AI, ...
            'Threshold',      P.DetThresh, ...
            'PsfFunPar',      {[0.1; 1.0; 1.5; 2.6]}, ...
            'RemoveEdgeDist', R + 2, ...
            'ReFind',         true);
    Cat   = AI.CatData;
    Bk    = AI.Back;  if isscalar(Bk); Bk = Bk + zeros(size(Im)); end
    Sub   = Im - Bk;
    clear AI Bk

    C = struct('NDet', 0, 'NSel', 0, 'SigSeed', NaN, 'S', emptyS(), ...
               'EAcc', 0, 'ECnt', 0, 'ENst', 0);
    if isempty(Cat.Catalog); return; end
    XY    = Cat.getCol({'X1', 'Y1'});
    XYpk  = Cat.getCol({'XPEAK', 'YPEAK'});
    SN    = Cat.getCol('SN_4');                   % S/N in the largest aperture
    AperCols = Cat.ColNames(startsWith(Cat.ColNames, 'FLUX_APER'));
    FluxA = Cat.getCol(AperCols{end});            % largest aperture
    M2    = Cat.getCol({'X2', 'Y2'});
    bad   = ~isfinite(XY(:,1)) | ~isfinite(XY(:,2));
    XY(bad,:) = XYpk(bad,:);
    if P.UseFlags && any(strcmp(Cat.ColNames, 'FLAGS'))
        Flags = Cat.getCol('FLAGS');
    else
        Flags = zeros(size(SN));
    end
    C.NDet = nnz(inTile(XY, x1, y1, xe(i), xe(i+1), ye(j), ye(j+1)));

    % candidates: owned by this tile, bright, clean, away from crop edges
    Sel = SN >= P.MinSN & Flags == 0 & FluxA > 0 & ...
          inTile(XY, x1, y1, xe(i), xe(i+1), ye(j), ye(j+1)) & ...
          XY(:,1) > R+2 & XY(:,1) < nx-R-1 & XY(:,2) > R+2 & XY(:,2) < ny-R-1;
    Cand = find(Sel);
    % isolation (neighbours searched in the whole crop, including the margin)
    Iso = true(size(Cand));
    for q = 1:numel(Cand)
        kk = Cand(q);
        dd = (XY(:,1)-XY(kk,1)).^2 + (XY(:,2)-XY(kk,2)).^2;
        Iso(q) = ~any(dd < P.IsoRad^2 & dd > 0 & FluxA > P.IsoFluxFrac*FluxA(kk));
    end
    Cand = Cand(Iso);
    % saturation / non-linearity
    Pk = zeros(size(Cand));
    for q = 1:numel(Cand)
        xi = round(XY(Cand(q),1)); yi = round(XY(Cand(q),2));
        Pk(q) = max(Im(yi-2:yi+2, xi-2:xi+2), [], 'all');
    end
    Cand = Cand(Pk < P.MaxPeakFrac*Satur);
    C.NSel = numel(Cand);
    C.SigSeed = sqrt(median(mean(M2(Cand,:), 2), 'omitnan'));
    if isempty(Sig) || isempty(Cand); return; end

    % moments (crop coordinates), then shift to full-image coordinates
    S = measureAll(Sub, XY(Cand,:), 'StampRad', R, ...
                   'SigWin', Sig*P.WinFactor, 'SigCore', Sig*P.CoreFactor);
    keep = isfinite(S.sigma) & isfinite(S.xcore) & isfinite(S.xbroad);
    fn = fieldnames(S);
    for f = 1:numel(fn); S.(fn{f}) = S.(fn{f})(keep); end
    S.sn   = SN(Cand(keep));
    S.crop = repmat((j-1)*numel(xe(1:end-1)) + i, nnz(keep), 1);
    xcC = S.xcore;  ycC = S.ycore;                % crop coords for ePSF
    S.xcore  = S.xcore  + x1 - 1;  S.ycore  = S.ycore  + y1 - 1;
    S.xbroad = S.xbroad + x1 - 1;  S.ybroad = S.ybroad + y1 - 1;
    C.S = S;

    % ePSF accumulation (stars passing the global size cut, top S/N per cell)
    if ~isempty(Psf)
        nPx = numel(Psf.pex) - 1;  nPy = numel(Psf.pey) - 1;
        nG  = 2*Psf.Hs*Psf.Ov + 1;
        C.EAcc = zeros(nG, nG, nPx*nPy);  C.ECnt = C.EAcc;  C.ENst = zeros(nPx*nPy, 1);
        good = S.sigma > 0.5*Sig & S.sigma < 3*Sig & S.flux > 0;
        cx = discretize(S.xcore, Psf.pex);  cy = discretize(S.ycore, Psf.pey);
        cc = (cy-1)*nPx + cx;
        quota = ceil(P.MaxStarsPsf / max(1, prod(P.CropGrid)/(nPx*nPy)));  % per crop
        for c = unique(cc(good & isfinite(cc)))'
            in = find(good & cc == c);
            [~, o] = sort(S.sn(in), 'descend');  in = in(o(1:min(end, quota)));
            [A, N] = accumEPSF(Sub, xcC(in), ycC(in), S.flux(in), 'HalfSize', Psf.Hs, 'Oversamp', Psf.Ov);
            C.EAcc(:,:,c) = A;  C.ECnt(:,:,c) = N;  C.ENst(c) = numel(in);
        end
    end
end

function f = inTile(XY, x1, y1, xa, xb, ya, yb)
% true if crop-coordinate XY lies in the interior tile (xa, xb] x (ya, yb] of the full image
    X = XY(:,1) + x1 - 1;  Y = XY(:,2) + y1 - 1;
    f = X > xa + 0.5 & X <= xb + 0.5 & Y > ya + 0.5 & Y <= yb + 0.5;
end

function S = emptyS()
    e = zeros(0,1);
    S = struct('xcore',e,'ycore',e,'xbroad',e,'ybroad',e,'sigma',e,'e1',e,'e2',e, ...
               'c1',e,'c2',e,'t1',e,'t2',e,'flux',e,'sn',e,'crop',e);
end

function S = catStruct(Sc)
% Concatenate a cell array of per-crop star structs.
    S = emptyS();  fn = fieldnames(S);
    for f = 1:numel(fn)
        v = cellfun(@(s) s.(fn{f})(:), Sc, 'UniformOutput', false);
        S.(fn{f}) = vertcat(v{:});
    end
end

function [Acc, Cnt] = accumEPSF(Sub, X, Y, Flux, Args)
    % Accumulate oversampled effective-PSF (ePSF) sums from a set of stars.
    %   Each pixel of each star, normalized by the star's flux, is dropped
    %   on a fine grid (Oversamp sub-pixels per pixel) at its offset from
    %   the star centroid. Sums and sample counts are returned separately so
    %   results from several crops can be added before normalizing:
    %   ePSF = Acc ./ max(Cnt, 1).
    % Input  : - Background-subtracted image (2D, double). Rows = Y, columns = X.
    %          - Vector of star X centroids [pix], in the image coordinates.
    %          - Vector of star Y centroids [pix].
    %          - Vector of star fluxes [ADU] used for normalization.
    %          * ...,key,val,...
    %            'HalfSize' - ePSF half-size [pix]. The fine grid covers
    %                   -HalfSize..+HalfSize pixels. Default is 6.
    %            'Oversamp' - Oversampling factor (sub-pixels per pixel).
    %                   Default is 4.
    % Output : - Matrix (nG x nG, nG = 2*HalfSize*Oversamp+1) of summed
    %            flux-normalized pixel values. Rows = Y offset, columns = X offset.
    %          - Matrix (nG x nG) of the number of samples in each fine pixel.
    % Example: [Acc, Cnt] = accumEPSF(Sub, S.xcore, S.ycore, S.flux, 'HalfSize', 6, 'Oversamp', 4);
    %          E = Acc ./ max(Cnt, 1);
    arguments
        Sub
        X(:,1)
        Y(:,1)
        Flux(:,1)
        Args.HalfSize(1,1) = 6;
        Args.Oversamp(1,1) = 4;
    end
    Hs = Args.HalfSize;  Ov = Args.Oversamp;  xc = X;  yc = Y;  flux = Flux;
    nG = 2*Hs*Ov + 1;  Acc = zeros(nG);  Cnt = zeros(nG);
    R = Hs + 1;
    for q = 1:numel(xc)
        xi = round(xc(q)); yi = round(yc(q));
        if xi-R < 1 || yi-R < 1 || xi+R > size(Sub,2) || yi+R > size(Sub,1); continue; end
        [X, Y] = meshgrid(xi-R:xi+R, yi-R:yi+R);
        V  = Sub(yi-R:yi+R, xi-R:xi+R) / flux(q);
        jx = round((X - xc(q))*Ov) + Hs*Ov + 1;  jy = round((Y - yc(q))*Ov) + Hs*Ov + 1;
        k  = jx >= 1 & jx <= nG & jy >= 1 & jy <= nG;
        idx = sub2ind([nG nG], jy(k), jx(k));
        Acc(:) = Acc(:) + accumarray(idx, V(k), [nG*nG 1]);
        Cnt(:) = Cnt(:) + accumarray(idx, 1,    [nG*nG 1]);
    end
end