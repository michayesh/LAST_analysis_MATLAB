function [Result, Aux] = measure_ACF_fwhm_fits(Input, Args)
%MEASURE_ACF_FWHM_FITS  Measure the PSF FWHM of a star-field FITS image with the
%   auto-correlation-function (ACF) method (imUtil.psf.fwhm_fromACF), over
%   the whole image or over sub-regions (CCDSEC), after removing saturated
%   and blended stars.
%
%   The ACF method works on pixels, not on a star list, so stars are
%   "filtered out" by: (1) detecting all sources, (2) flagging saturated
%   and blended ones, and (3) replacing their pixels (plus a safety margin)
%   with local background + Gaussian noise before the ACF is computed.
%   Isolated, unsaturated stars are left untouched.
%
%   Syntax
%     Result        = measure_ACF_fwhm_fits(FileName)
%     [Result, Aux] = measure_ACF_fwhm_fits(FileName, Name=Value, ...)
%     [Result, Aux] = measure_ACF_fwhm_fits(ImageMatrix, Name=Value, ...)
%
%   Region selection (default: whole image)
%     CCDSEC   N-by-4 matrix, one region per row: [Xmin Xmax Ymin Ymax]
%              (X = column, Y = row), same convention as fwhm_fromACF.
%     Grid     [Ntx Nty] - split the image into Ntx-by-Nty equal tiles
%              (scalar N means N-by-N).
%
%   FITS input
%     HDU      HDU number to read (1 = primary). [] = first HDU that holds
%              a 2-D image (works for fpack/tile-compressed files too).
%
%   Source detection
%     BackBlock     Block size [pix] of the background/noise map. Default 64.
%     SmoothSigma   Gaussian smoothing sigma [pix] for detection. Default 1.5.
%     DetectNsigma  Detection threshold in sigma. Default 5.
%     MinArea       Minimum detected area [pix]. Default 5.
%
%   Saturation
%     SatLevel      Saturation level [ADU]. [] = FITS keyword SATURATE
%                   (or SATLEVEL), otherwise 60000. Default [].
%     NonLinFrac    Stars whose peak >= NonLinFrac*SatLevel are rejected
%                   (removes the non-linear regime too). Default 0.9.
%     SatMaskGrow   Mask dilation radius around saturated stars [pix].
%                   Default 15.
%
%   Blending
%     MinSep          Sources whose nearest neighbour is closer than this
%                     [pix] are rejected as blended. Default 10.
%     PeakProminence  A source with more than one peak of at least this
%                     prominence [sigma] is rejected as blended. Default 3.
%     MaskGrow        Mask dilation radius around blended stars [pix].
%                     Default 5.
%
%   ACF / FWHM (passed to fwhm_fromACF)
%     CorrFrac        Default 0.65.
%     Nsigma0         Default 10.
%     MaxRadius       Default 50 (clipped to fit each region).
%     Step            Default 1.
%     MedianSatLevel  fwhm_fromACF 'SatLevel' (median-based check). Default 30000.
%     FlattenBack     Subtract the background map (keeping its median level)
%                     so gradients do not leak into the ACF. Default true.
%     Engine          "auto" | "astropack" | "local". "auto" uses AstroPack's
%                     imUtil.psf.fwhm_fromACF if it (and its dependencies) are
%                     on the path, otherwise a self-contained re-implementation.
%
%   Misc
%     MinStars   Regions with fewer good stars get LowStars=true. Default 3.
%     PixScale   [arcsec/pix]; if given, adds FWHM_arcsec column. Default [].
%     RandSeed   Seed for the noise used to fill masked pixels. Default 1.
%     Plot       Show diagnostic figures. Default true.
%     Verbose    Print a summary. Default true.
%
%   Output
%     Result  table, one row per region: Region, Xmin, Xmax, Ymin, Ymax,
%             Xcenter, Ycenter, FWHM, Status, NgoodStars, NsatStars,
%             NblendStars, MaskedFrac, LowStars [, FWHM_arcsec].
%             Note: as in fwhm_fromACF, FWHM is an ACF-based quantity that
%             scales like the FWHM, not a formal Gaussian FWHM.
%     Aux     struct with Sources (table), Mask, Back, Std, CleanImage,
%             SatLevel, Engine, Sections.
%
%   Examples
%     R = measure_ACF_fwhm_fits('field.fits');
%     R = measure_ACF_fwhm_fits('field.fits', Grid=[3 3]);
%     R = measure_ACF_fwhm_fits('field.fits', CCDSEC=[1 2000 1 2000; 2001 4000 1 2000]);
%     [R, A] = measure_ACF_fwhm_fits('field.fits', Grid=4, SatLevel=50000, MinSep=15);
%
%   Requires: Image Processing Toolbox, Statistics and Machine Learning Toolbox.

arguments
    Input
    % --- input / regions
    Args.HDU                = [];
    Args.CCDSEC             = [];
    Args.Grid               = [];
    % --- detection
    Args.BackBlock      (1,1) double {mustBePositive}    = 64;
    Args.SmoothSigma    (1,1) double {mustBePositive}    = 1.5;
    Args.DetectNsigma   (1,1) double {mustBePositive}    = 5;
    Args.MinArea        (1,1) double {mustBeNonnegative} = 5;
    % --- saturation
    Args.SatLevel                                        = [];
    Args.NonLinFrac     (1,1) double {mustBePositive}    = 0.9;
    Args.SatMaskGrow    (1,1) double {mustBeNonnegative} = 15;
    % --- blending
    Args.MinSep         (1,1) double {mustBeNonnegative} = 10;
    Args.PeakProminence (1,1) double {mustBePositive}    = 3;
    Args.MaskGrow       (1,1) double {mustBeNonnegative} = 5;
    % --- ACF
    Args.CorrFrac       (1,1) double = 0.65;
    Args.Nsigma0        (1,1) double = 10;
    Args.MaxRadius      (1,1) double {mustBePositive} = 50;
    Args.Step           (1,1) double {mustBePositive} = 1;
    Args.MedianSatLevel (1,1) double = 30000;
    Args.FlattenBack    (1,1) logical = true;
    Args.Engine         (1,1) string {mustBeMember(Args.Engine,["auto","astropack","local"])} = "auto";
    % --- misc
    Args.MinStars       (1,1) double = 3;
    Args.PixScale                    = [];
    Args.RandSeed       (1,1) double = 1;
    Args.Plot           (1,1) logical = true;
    Args.Verbose        (1,1) logical = true;
end

%% 1. Load image -----------------------------------------------------------
if isnumeric(Input) || islogical(Input)
    Image  = double(Input);
    HdrSat = NaN;
    Name   = "input matrix";
else
    [Image, HdrSat] = local_read_fits(char(Input), Args.HDU);
    Name = string(Input);
end
if ~ismatrix(Image)
    warning('measure_ACF_fwhm_fits:Cube', 'Input is a cube; using the first plane.');
    Image = Image(:,:,1);
end
[Ny, Nx] = size(Image);

%% 2. Saturation level -----------------------------------------------------
if ~isempty(Args.SatLevel)
    SatLevel = Args.SatLevel;
elseif isfinite(HdrSat)
    SatLevel = HdrSat;
else
    SatLevel = 60000;
    warning('measure_ACF_fwhm_fits:NoSatLevel', ...
        'No SATURATE keyword found; using SatLevel = %g. Set SatLevel explicitly if wrong.', SatLevel);
end

%% 3. Background and noise maps -------------------------------------------
[Back, Std] = local_background(Image, Args.BackBlock);

%% 4. Source detection -----------------------------------------------------
Sub = Image - Back;
Sub(~isfinite(Sub)) = 0;
Smooth    = imgaussfilt(Sub, Args.SmoothSigma);
SmoothStd = Std ./ (2*sqrt(pi)*Args.SmoothSigma);   % noise after Gaussian smoothing

Detect = Smooth > Args.DetectNsigma .* SmoothStd;
Detect = bwareaopen(Detect, max(1, round(Args.MinArea)), 8);
CC     = bwconncomp(Detect, 8);
Lab    = labelmatrix(CC);
Nsrc   = CC.NumObjects;

PropsW = regionprops(CC, Sub,   'WeightedCentroid');
PropsI = regionprops(CC, Image, 'Area', 'MaxIntensity');
XY     = reshape([PropsW.WeightedCentroid], 2, []).';
Area   = reshape([PropsI.Area], [], 1);
Peak   = reshape([PropsI.MaxIntensity], [], 1);

% number of significant peaks inside each detected component
H        = Args.PeakProminence * median(SmoothStd(:), 'omitnan');
PeakMap  = imextendedmax(Smooth, H, 8) & Detect;
PCC      = bwconncomp(PeakMap, 8);
PeakLab  = cellfun(@(p) double(Lab(p(1))), PCC.PixelIdxList);
PeakLab  = PeakLab(PeakLab > 0);
Npeaks   = accumarray(PeakLab(:), 1, [Nsrc 1]);

% nearest-neighbour distance
if Nsrc >= 2
    [~, D] = knnsearch(XY, XY, 'K', 2);
    NNdist = D(:,2);
else
    NNdist = inf(Nsrc, 1);
end

%% 5. Flag saturated and blended stars ------------------------------------
IsSat   = Peak >= Args.NonLinFrac * SatLevel;
IsBlend = (Npeaks > 1) | (NNdist < Args.MinSep);
IsGood  = ~IsSat & ~IsBlend;

%% 6. Mask bad stars and fill with background + noise ---------------------
MaskSat = ismember(Lab, find(IsSat));
if Args.SatMaskGrow > 0 && any(MaskSat(:))
    MaskSat = imdilate(MaskSat, strel('disk', round(Args.SatMaskGrow), 0));
end
MaskBl = ismember(Lab, find(IsBlend & ~IsSat));
if Args.MaskGrow > 0 && any(MaskBl(:))
    MaskBl = imdilate(MaskBl, strel('disk', round(Args.MaskGrow), 0));
end
Mask = MaskSat | MaskBl;

Clean = Image;
Nmask = nnz(Mask);
if Nmask > 0
    RS = RandStream('mt19937ar', 'Seed', Args.RandSeed);
    Clean(Mask) = Back(Mask) + Std(Mask) .* randn(RS, Nmask, 1);
end

if Args.FlattenBack
    ImForACF = Clean - Back + median(Back(:), 'omitnan');
else
    ImForACF = Clean;
end

%% 7. Choose FWHM engine ---------------------------------------------------
HasAP = ~isempty(which('imUtil.psf.fwhm_fromACF'))   && ...
        ~isempty(which('tools.math.stat.std_mad'))    && ...
        ~isempty(which('imUtil.psf.radialProfile'))   && ...
        ~isempty(which('tools.interp.interp1crossVal'));
switch Args.Engine
    case "auto"
        if HasAP, Engine = "astropack"; else, Engine = "local"; end
    case "astropack"
        if ~HasAP
            error('measure_ACF_fwhm_fits:NoAstroPack', ...
                'AstroPack functions (imUtil.psf.fwhm_fromACF and dependencies) not found on the path.');
        end
        Engine = "astropack";
    otherwise
        Engine = "local";
end
if Engine == "astropack"
    FwhmFun = @(varargin) imUtil.psf.fwhm_fromACF(varargin{:});
else
    FwhmFun = @local_fwhm_fromACF;
end

%% 8. Regions and FWHM measurement -----------------------------------------
[Sec, GridSize] = local_make_sections(Nx, Ny, Args.CCDSEC, Args.Grid);
Nsec = size(Sec, 1);

[FWHM, Xc, Yc, MaskedFrac] = deal(nan(Nsec, 1));
[Ngood, Nsat, Nblend]      = deal(zeros(Nsec, 1));
[Status, LowStars]         = deal(false(Nsec, 1));

for k = 1:Nsec
    S = Sec(k,:);
    Xc(k) = (S(1) + S(2)) / 2;
    Yc(k) = (S(3) + S(4)) / 2;

    InSec = XY(:,1) >= S(1) & XY(:,1) <= S(2) & XY(:,2) >= S(3) & XY(:,2) <= S(4);
    Ngood(k)  = nnz(InSec & IsGood);
    Nsat(k)   = nnz(InSec & IsSat);
    Nblend(k) = nnz(InSec & IsBlend & ~IsSat);
    MaskedFrac(k) = mean(Mask(S(3):S(4), S(1):S(2)), 'all');
    LowStars(k)   = Ngood(k) < Args.MinStars;

    RegSize = [S(4)-S(3)+1, S(2)-S(1)+1];
    MaxRad  = min(Args.MaxRadius, floor(min(RegSize)/2) - 1);
    if MaxRad < 3*Args.Step
        warning('measure_ACF_fwhm_fits:SmallRegion', 'Region %d is too small; skipped.', k);
        continue
    end

    try
        [F, ~, Info] = FwhmFun(ImForACF, 'CCDSEC', S, ...
            'CorrFrac', Args.CorrFrac, 'Nsigma0', Args.Nsigma0, ...
            'MaxRadius', MaxRad, 'Step', Args.Step, ...
            'SatLevel', Args.MedianSatLevel);
        FWHM(k)   = F;
        Status(k) = logical(Info.Status) && isfinite(F);
    catch ME
        warning('measure_ACF_fwhm_fits:ACFfailed', 'Region %d: %s', k, ME.message);
    end
end

Result = table((1:Nsec).', Sec(:,1), Sec(:,2), Sec(:,3), Sec(:,4), Xc, Yc, ...
    FWHM, Status, Ngood, Nsat, Nblend, MaskedFrac, LowStars, ...
    'VariableNames', {'Region','Xmin','Xmax','Ymin','Ymax','Xcenter','Ycenter', ...
    'FWHM','Status','NgoodStars','NsatStars','NblendStars','MaskedFrac','LowStars'});
if ~isempty(Args.PixScale)
    Result.FWHM_arcsec = Result.FWHM .* Args.PixScale;
end

Sources = table(XY(:,1), XY(:,2), Area, Peak, Npeaks, NNdist, IsSat, IsBlend, IsGood, ...
    'VariableNames', {'X','Y','Area','Peak','Npeaks','NNdist','IsSat','IsBlend','IsGood'});

Aux = struct('Sources', Sources, 'Mask', Mask, 'Back', Back, 'Std', Std, ...
    'CleanImage', Clean, 'SatLevel', SatLevel, 'Engine', Engine, 'Sections', Sec);

%% 9. Report and plots ------------------------------------------------------
if Args.Verbose
    fprintf('\n%s  [%d x %d pix]\n', Name, Nx, Ny);
    fprintf('Engine: %s | SatLevel: %g | sources: %d (good %d, saturated %d, blended %d)\n', ...
        Engine, SatLevel, Nsrc, nnz(IsGood), nnz(IsSat), nnz(IsBlend & ~IsSat));
    fprintf('Masked pixels: %.2f%%\n\n', 100*Nmask/numel(Mask));
    disp(Result)
end

if Args.Plot
    local_plot(Image, Sec, Result, Sources, GridSize, Name);
end
end


%% ========================================================================
%  Local functions
%% ========================================================================

function [Image, SatHdr] = local_read_fits(FileName, HDU)
% Read a 2-D image HDU and look for a saturation keyword.
import matlab.io.*
if ~isfile(FileName)
    error('measure_ACF_fwhm_fits:FileNotFound', 'File not found: %s', FileName);
end
fptr    = fits.openFile(FileName);
Cleaner = onCleanup(@() fits.closeFile(fptr));

if isempty(HDU)
    HDU  = NaN;
    Nhdu = fits.getNumHDUs(fptr);
    for h = 1:Nhdu
        HduType = fits.movAbsHDU(fptr, h);
        if strcmp(HduType, 'IMAGE_HDU')
            Sz = fits.getImgSize(fptr);
            if numel(Sz) >= 2 && all(Sz > 0)
                HDU = h;
                break
            end
        end
    end
    if isnan(HDU)
        error('measure_ACF_fwhm_fits:NoImage', 'No 2-D image HDU found in %s', FileName);
    end
end

fits.movAbsHDU(fptr, HDU);
Image = double(fits.readImg(fptr));   % rows = NAXIS2 (Y), columns = NAXIS1 (X)

SatHdr = NaN;
for Key = ["SATURATE", "SATLEVEL", "SATURATN"]
    try
        V = str2double(fits.readKey(fptr, char(Key)));
        if isfinite(V)
            SatHdr = V;
            break
        end
    catch
        % keyword not present
    end
end
end


function [Back, Std] = local_background(Image, Block)
% Block-wise sigma-clipped median background and robust (MAD) noise maps.
[Ny, Nx] = size(Image);
Nby = max(1, round(Ny / Block));
Nbx = max(1, round(Nx / Block));
Ye  = round(linspace(0, Ny, Nby + 1));
Xe  = round(linspace(0, Nx, Nbx + 1));

B = nan(Nby, Nbx);
S = nan(Nby, Nbx);
for i = 1:Nby
    for j = 1:Nbx
        T = Image(Ye(i)+1:Ye(i+1), Xe(j)+1:Xe(j+1));
        T = T(isfinite(T));
        if numel(T) < 10, continue; end
        for Iter = 1:2                         % two passes of 3-sigma clipping
            M = median(T);
            Sd = 1.4826 * median(abs(T - M));
            if Sd > 0
                T = T(abs(T - M) < 3*Sd);
            end
        end
        B(i,j) = median(T);
        S(i,j) = 1.4826 * median(abs(T - B(i,j)));
    end
end

% fill empty cells, reject star-dominated cells, and fix zero noise
B = local_fill(B);
S = local_fill(S);
if Nby >= 3 && Nbx >= 3
    B = medfilt2(B, [3 3], 'symmetric');
    S = medfilt2(S, [3 3], 'symmetric');
end
Pos = S(S > 0);
if isempty(Pos), Pos = 1; end
S(~(S > 0)) = median(Pos);

Back = imresize(B, [Ny Nx], 'bilinear');
Std  = imresize(S, [Ny Nx], 'bilinear');
end


function A = local_fill(A)
if any(isnan(A(:)))
    A = fillmissing(A, 'nearest', 1);
    A = fillmissing(A, 'nearest', 2);
    A(isnan(A)) = median(A(:), 'omitnan');
end
end


function [Sec, GridSize] = local_make_sections(Nx, Ny, CCDSEC, Grid)
% Build the list of regions [Xmin Xmax Ymin Ymax].
GridSize = [];
if ~isempty(CCDSEC) && ~isempty(Grid)
    error('measure_ACF_fwhm_fits:RegionArgs', 'Specify either CCDSEC or Grid, not both.');
end

if ~isempty(CCDSEC)
    if size(CCDSEC, 2) ~= 4
        error('measure_ACF_fwhm_fits:CCDSEC', 'CCDSEC must be N-by-4: [Xmin Xmax Ymin Ymax].');
    end
    Sec = round(double(CCDSEC));
    Orig = Sec;
    Sec(:,[1 2]) = min(max(Sec(:,[1 2]), 1), Nx);
    Sec(:,[3 4]) = min(max(Sec(:,[3 4]), 1), Ny);
    if ~isequal(Sec, Orig)
        warning('measure_ACF_fwhm_fits:CCDSECclipped', 'Some CCDSEC rows were clipped to the image bounds.');
    end
    if any(Sec(:,2) <= Sec(:,1) | Sec(:,4) <= Sec(:,3))
        error('measure_ACF_fwhm_fits:CCDSEC', 'Each CCDSEC row needs Xmin<Xmax and Ymin<Ymax.');
    end

elseif ~isempty(Grid)
    if isscalar(Grid), Grid = [Grid Grid]; end
    Grid = round(Grid(:).');
    Xe = round(linspace(0, Nx, Grid(1) + 1));
    Ye = round(linspace(0, Ny, Grid(2) + 1));
    Sec = zeros(prod(Grid), 4);
    k = 0;
    for j = 1:Grid(2)          % Y tiles (outer)
        for i = 1:Grid(1)      % X tiles (inner, fastest)
            k = k + 1;
            Sec(k,:) = [Xe(i)+1, Xe(i+1), Ye(j)+1, Ye(j+1)];
        end
    end
    GridSize = Grid;

else
    Sec = [1 Nx 1 Ny];
end
end


function [FWHM, Nstars, Info, ACF] = local_fwhm_fromACF(Image, Args)
% Self-contained re-implementation of imUtil.psf.fwhm_fromACF (used when
% AstroPack is not on the path). Same algorithm; radial binning may differ
% slightly from imUtil.psf.radialProfile, so absolute values can differ a
% little from the AstroPack version.
arguments
    Image
    Args.CCDSEC    = [];
    Args.CorrFrac  = 0.65;
    Args.Nsigma0   = 10;
    Args.MaxRadius = 200;
    Args.Step      = 1;
    Args.SatLevel  = 30000;
end
Nstars = NaN;
FWHM   = NaN;
ACF    = [];
Info   = struct('Status', false, 'Rad', [], 'CumVal', []);

if ~isempty(Args.CCDSEC)
    S = Args.CCDSEC(1,:);
    Image = Image(S(3):S(4), S(1):S(2));
end
Image = single(Image);

Med = median(Image, 'all', 'omitnan');
if Med > Args.SatLevel
    return
end

Image(isnan(Image)) = 0;
Image = Image - Med;
Std   = 1.4826 * median(abs(Image(:) - median(Image(:))));
Image(Image < Args.Nsigma0 .* Std) = 0;

F   = fft2(Image);
ACF = fftshift(real(ifft2(F .* conj(F))));

[Ny, Nx] = size(ACF);
Cy = floor(Ny/2) + 1;                    % zero-lag pixel after fftshift
Cx = floor(Nx/2) + 1;
MaxR = min([Args.MaxRadius, Cy-1, Cx-1, Ny-Cy, Nx-Cx]);
MaxR = floor(MaxR);
if MaxR < 2*Args.Step
    return
end

Cut = double(ACF(Cy-MaxR:Cy+MaxR, Cx-MaxR:Cx+MaxR));
[X, Y] = meshgrid(-MaxR:MaxR);
R = sqrt(X.^2 + Y.^2);
Edges = 0:Args.Step:MaxR;
Bin = discretize(R, Edges);
Ok  = ~isnan(Bin);
Nb  = numel(Edges) - 1;
MeanV = accumarray(Bin(Ok), Cut(Ok), [Nb 1], @mean, NaN);
MeanR = accumarray(Bin(Ok), R(Ok),   [Nb 1], @mean, NaN);

% skip the zero-lag bin (noise spike), as the original does
Rad   = MeanR(2:end);
MeanV = MeanV(2:end);
Keep  = ~isnan(Rad);
Rad   = Rad(Keep);
MeanV = MeanV(Keep);
MeanV(isnan(MeanV)) = 0;

CumVal = cumsum(MeanV);
CumVal = CumVal ./ CumVal(end);
if isempty(CumVal) || any(~isfinite(CumVal))
    return
end
CumVal = CumVal(:) + (1:numel(CumVal)).' .* 1e-7;   % enforce monotonic

Idx = find(CumVal >= Args.CorrFrac, 1, 'first');
if isempty(Idx)
    return
elseif Idx == 1
    FWHM = Rad(1) * Args.CorrFrac / CumVal(1);
else
    FWHM = interp1(CumVal(Idx-1:Idx), Rad(Idx-1:Idx), Args.CorrFrac);
end

Info.Rad    = Rad;
Info.CumVal = CumVal;
Info.Status = isfinite(FWHM);
end


function local_plot(Image, Sec, Result, Sources, GridSize, Name)
% Image with sources and regions, plus a FWHM map for multi-region runs.
V = Image(isfinite(Image));
if numel(V) > 1e6
    V = V(randperm(numel(V), 1e6));
end
M = median(V);
S = 1.4826 * median(abs(V - M));
if ~(S > 0), S = std(V); end
if ~(S > 0), S = 1; end

figure('Name', 'FWHM: image and sources', 'Color', 'w');
imagesc(Image, [M - 2*S, M + 15*S]);
axis image; set(gca, 'YDir', 'normal'); colormap(gca, gray); colorbar;
hold on
G = Sources.IsGood;
Sa = Sources.IsSat;
Bl = Sources.IsBlend & ~Sources.IsSat;
plot(Sources.X(G),  Sources.Y(G),  'o', 'Color', [0 0.8 0], 'MarkerSize', 6, 'DisplayName', 'good');
plot(Sources.X(Sa), Sources.Y(Sa), 'o', 'Color', [1 0 0],   'MarkerSize', 12, 'LineWidth', 1.2, 'DisplayName', 'saturated');
plot(Sources.X(Bl), Sources.Y(Bl), 's', 'Color', [1 0.6 0], 'MarkerSize', 9,  'LineWidth', 1.2, 'DisplayName', 'blended');
for k = 1:size(Sec, 1)
    Sk = Sec(k,:);
    rectangle('Position', [Sk(1)-0.5, Sk(3)-0.5, Sk(2)-Sk(1)+1, Sk(4)-Sk(3)+1], 'EdgeColor', 'c');
    text(Result.Xcenter(k), Result.Ycenter(k), sprintf('%.2f', Result.FWHM(k)), ...
        'Color', 'y', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
legend('Location', 'northeastoutside');
xlabel('X [pix]'); ylabel('Y [pix]');
title(sprintf('%s  (ACF FWHM per region)', Name), 'Interpreter', 'none');
hold off

if height(Result) > 1
    figure('Name', 'FWHM map', 'Color', 'w');
    if ~isempty(GridSize)
        Map = reshape(Result.FWHM, GridSize(1), GridSize(2)).';
        Xv  = Result.Xcenter(1:GridSize(1));
        Yv  = Result.Ycenter(1:GridSize(1):end);
        imagesc(Xv, Yv, Map);
        axis image; set(gca, 'YDir', 'normal'); colorbar;
        for k = 1:height(Result)
            text(Result.Xcenter(k), Result.Ycenter(k), ...
                sprintf('%.2f\n(%d)', Result.FWHM(k), Result.NgoodStars(k)), ...
                'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
        end
    else
        scatter(Result.Xcenter, Result.Ycenter, 200, Result.FWHM, 'filled', 'MarkerEdgeColor', 'k');
        axis equal; colorbar; grid on;
    end
    xlabel('X [pix]'); ylabel('Y [pix]');
    title('ACF FWHM (number of good stars)');
end
end