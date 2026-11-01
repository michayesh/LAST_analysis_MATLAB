function [Result, Aux] = measure_wholeImageFWHM_from_fits(Input, Args)
%MEASURE_WHOLEIMAGEFWHM_FROM_FITS  Measure the PSF FWHM of a star-field FITS
%   image with the AstroPack whole-image estimators (imUtil.psf.fwhm_fromACF,
%   fwhm_fromMoments, fwhm_fromBank), over the whole image or over
%   sub-regions (CCDSEC), after removing saturated and blended stars.
%
%   The estimators work on pixels, so stars are "filtered out" by:
%   (1) detecting all sources, (2) flagging saturated and blended ones, and
%   (3) replacing their pixels (plus a safety margin) with local background
%   + Gaussian noise before the FWHM is measured. Isolated, unsaturated
%   stars are left untouched. All methods run on the same cleaned image and
%   the same regions, so their results can be compared directly.
%
%   Syntax
%     Result        = measure_wholeImageFWHM_from_fits(FileName)
%     [Result, Aux] = measure_wholeImageFWHM_from_fits(FileName, Name=Value, ...)
%     [Result, Aux] = measure_wholeImageFWHM_from_fits(ImageMatrix, Name=Value, ...)
%
%   Method
%     Method   "acf" | "moments" | "bank" | "all", or a string array of
%              several, e.g. ["acf","moments"]. Default "acf".
%                "acf"     - imUtil.psf.fwhm_fromACF: auto-correlation;
%                            best for large / out-of-focus PSFs. The value
%                            scales like the FWHM but is not a formal FWHM.
%                "moments" - imUtil.psf.fwhm_fromMoments: median 2nd
%                            moments of detected stars; best near focus.
%                "bank"    - imUtil.psf.fwhm_fromBank: Gaussian template
%                            bank, most detections above S/N.
%                "all"     - all three.
%              The first method listed is the "primary" one (FWHM and
%              Status columns, plots). "moments" and "bank" need AstroPack;
%              "acf" falls back to a built-in copy if AstroPack is missing.
%
%   Region selection (default: whole image)
%     CCDSEC   N-by-4 matrix, one region per row: [Xmin Xmax Ymin Ymax]
%              (X = column, Y = row), same convention as the imUtil.psf
%              functions.
%     Grid     [Ntx Nty] - split the image into Ntx-by-Nty equal tiles
%              (scalar N means N-by-N).
%
%   FITS input
%     HDU      HDU number to read (1 = primary). [] = first HDU that holds
%              a 2-D image (works for fpack/tile-compressed files too).
%
%   Source detection (for masking)
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
%   Common
%     FlattenBack     Subtract the background map (keeping its median level)
%                     so gradients do not bias the estimators. Default true.
%
%   ACF options (fwhm_fromACF)
%     CorrFrac        Default 0.65.
%     Nsigma0         Default 10.
%     MaxRadius       Default 50 (clipped to fit each region).
%     Step            Default 1.
%     MedianSatLevel  fwhm_fromACF 'SatLevel' (median-based check). Default 30000.
%     ACFEngine       "auto" | "astropack" | "local". "auto" uses AstroPack's
%                     imUtil.psf.fwhm_fromACF if it (and its dependencies) are
%                     on the path, otherwise a self-contained re-implementation.
%
%   Moments options (fwhm_fromMoments)
%     MomMinSN        Detection S/N threshold. Default 30.
%     MomRadius       Radius for the moment calculation [pix]
%                     (recommended ~1.7 FWHM). Default 12.
%
%   Bank options (fwhm_fromBank)
%     BankMinSN       Detection S/N threshold. Default 50.
%     BankMethod      "bisec" | "maxndet" | "maxndetinterp". Default "bisec".
%     BankSigmaVec    Gaussian sigmas of the bank [pix]. Default logspace(0,2,5).'.
%     BankMinStars    Minimum stars for an estimate. Default 5.
%     BankMaxIter     Iterations of the bisec method. Default 6.
%     BankUseBackMap  Pass this function's background and variance maps to
%                     fwhm_fromBank instead of its internal estimate (which
%                     assumes Poisson noise with gain 1). Default true.
%                     fwhm_fromBank is always called with PixScale=1, so
%                     its FWHM is in pixels here.
%
%   Misc
%     MinStars   Regions with fewer good stars get LowStars=true. Default 3.
%     PixScale   [arcsec/pix]; if given, adds *_arcsec columns. Default [].
%     RandSeed   Seed for the noise used to fill masked pixels. Default 1.
%     Plot       Show diagnostic figures. Default true.
%     MapFontSize  Font size of the numbers in the FWHM map cells. Default 7.
%     Verbose    Print the run stages (with elapsed time) and a summary
%                to the command window. Default true.
%
%   Output
%     Result  table, one row per region: Region, Xmin, Xmax, Ymin, Ymax,
%             Xcenter, Ycenter, FWHM, Status (primary method), NgoodStars,
%             NsatStars, NblendStars, MaskedFrac, LowStars, and for every
%             method M: FWHM_M, Status_M, Nused_M (stars used by the method;
%             NaN for acf) [, FWHM_arcsec, FWHM_M_arcsec]. FWHM in pixels.
%     Aux     struct with Sources (table), Mask, Back, Std, CleanImage,
%             SatLevel, Methods, ACFEngine, Sections.
%
%   Examples
%     R = measure_wholeImageFWHM_from_fits('field.fits');
%     R = measure_wholeImageFWHM_from_fits('field.fits', Method="moments", Grid=[3 3]);
%     R = measure_wholeImageFWHM_from_fits('field.fits', Method="all", Grid=4);
%     R = measure_wholeImageFWHM_from_fits('field.fits', Method=["bank","acf"], ...
%             CCDSEC=[1 2000 1 2000; 2001 4000 1 2000], PixScale=1.25);
%
%   Requires: Image Processing Toolbox, Statistics and Machine Learning
%   Toolbox; AstroPack for Method "moments" and "bank".

arguments
    Input
    % --- method
    Args.Method         (1,:) string {mustBeMember(Args.Method,["acf","moments","bank","all"])} = "acf";
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
    % --- common
    Args.FlattenBack    (1,1) logical = true;
    % --- ACF
    Args.CorrFrac       (1,1) double = 0.65;
    Args.Nsigma0        (1,1) double = 10;
    Args.MaxRadius      (1,1) double {mustBePositive} = 50;
    Args.Step           (1,1) double {mustBePositive} = 1;
    Args.MedianSatLevel (1,1) double = 30000;
    Args.ACFEngine      (1,1) string {mustBeMember(Args.ACFEngine,["auto","astropack","local"])} = "auto";
    % --- moments
    Args.MomMinSN       (1,1) double {mustBePositive} = 30;
    Args.MomRadius      (1,1) double {mustBePositive} = 12;
    % --- bank
    Args.BankMinSN      (1,1) double {mustBePositive} = 50;
    Args.BankMethod     (1,1) string {mustBeMember(Args.BankMethod,["bisec","maxndet","maxndetinterp"])} = "bisec";
    Args.BankSigmaVec   (:,1) double {mustBePositive} = logspace(0,2,5).';
    Args.BankMinStars   (1,1) double = 5;
    Args.BankMaxIter    (1,1) double {mustBePositive} = 6;
    Args.BankUseBackMap (1,1) logical = true;
    % --- misc
    Args.MinStars       (1,1) double = 3;
    Args.PixScale                    = [];
    Args.RandSeed       (1,1) double = 1;
    Args.Plot           (1,1) logical = true;
    Args.MapFontSize    (1,1) double {mustBePositive} = 7;
    Args.Verbose        (1,1) logical = true;
end

% stage messages on the command window (silenced with Verbose=false)
T0  = tic;
Log = @(varargin) local_log(Args.Verbose, T0, varargin{:});
Log('measure_wholeImageFWHM_from_fits: start');

%% 0. Methods and dependencies ---------------------------------------------
Log('[0/8] Checking methods and dependencies');
Methods = unique(lower(Args.Method), 'stable');
if any(Methods == "all")
    Methods = ["acf", "moments", "bank"];
end
Nm = numel(Methods);

has = @(names) all(cellfun(@(n) ~isempty(which(n)), names));
HasACF  = has({'imUtil.psf.fwhm_fromACF', 'tools.math.stat.std_mad', ...
               'imUtil.psf.radialProfile', 'tools.interp.interp1crossVal'});
HasMom  = has({'imUtil.psf.fwhm_fromMoments', 'imUtil.sources.findSources', ...
               'imUtil.image.moment2', 'tools.math.stat.rstd', 'imUtil.kernel2.gauss'});
HasBank = has({'imUtil.psf.fwhm_fromBank', 'imUtil.filter.filter2_snBank', ...
               'imUtil.image.local_maxima', 'imUtil.kernel2.gauss'});

if any(Methods == "moments") && ~HasMom
    error('measure_wholeImageFWHM_from_fits:NoAstroPack', ...
        'Method "moments" needs AstroPack (imUtil.psf.fwhm_fromMoments and dependencies) on the path.');
end
if any(Methods == "bank") && ~HasBank
    error('measure_wholeImageFWHM_from_fits:NoAstroPack', ...
        'Method "bank" needs AstroPack (imUtil.psf.fwhm_fromBank and dependencies) on the path.');
end
if Args.BankMethod == "maxndetinterp" && any(Methods == "bank") && isempty(which('tools.find.find_local_extremum'))
    error('measure_wholeImageFWHM_from_fits:NoAstroPack', ...
        'BankMethod "maxndetinterp" needs tools.find.find_local_extremum (AstroPack).');
end

switch Args.ACFEngine
    case "auto"
        if HasACF, ACFEngine = "astropack"; else, ACFEngine = "local"; end
    case "astropack"
        if ~HasACF && any(Methods == "acf")
            error('measure_wholeImageFWHM_from_fits:NoAstroPack', ...
                'AstroPack functions (imUtil.psf.fwhm_fromACF and dependencies) not found on the path.');
        end
        ACFEngine = "astropack";
    otherwise
        ACFEngine = "local";
end
if ACFEngine == "astropack"
    AcfFun = @(varargin) imUtil.psf.fwhm_fromACF(varargin{:});
else
    AcfFun = @local_fwhm_fromACF;
end
Log('      methods: %s (primary: %s) | ACF engine: %s', strjoin(Methods, ', '), Methods(1), ACFEngine);

%% 1. Load image -----------------------------------------------------------
Log('[1/8] Loading image');
if isnumeric(Input) || islogical(Input)
    Image  = double(Input);
    HdrSat = NaN;
    Name   = "input matrix";
else
    [Image, HdrSat] = local_read_fits(char(Input), Args.HDU);
    Name = string(Input);
end
if ~ismatrix(Image)
    warning('measure_wholeImageFWHM_from_fits:Cube', 'Input is a cube; using the first plane.');
    Image = Image(:,:,1);
end
[Ny, Nx] = size(Image);
Log('      %s  [%d x %d pix]', Name, Nx, Ny);

%% 2. Saturation level -----------------------------------------------------
Log('[2/8] Setting saturation level');
if ~isempty(Args.SatLevel)
    SatLevel = Args.SatLevel;
    SatSrc   = 'SatLevel argument';
elseif isfinite(HdrSat)
    SatLevel = HdrSat;
    SatSrc   = 'FITS header';
else
    SatLevel = 60000;
    SatSrc   = 'default';
    warning('measure_wholeImageFWHM_from_fits:NoSatLevel', ...
        'No SATURATE keyword found; using SatLevel = %g. Set SatLevel explicitly if wrong.', SatLevel);
end
Log('      SatLevel = %g (%s); stars with peak >= %g are rejected', ...
    SatLevel, SatSrc, Args.NonLinFrac*SatLevel);

%% 3. Background and noise maps -------------------------------------------
Log('[3/8] Computing background and noise maps (block %g pix)', Args.BackBlock);
[Back, Std] = local_background(Image, Args.BackBlock);
Log('      median background = %.1f, median noise = %.2f', ...
    median(Back(:), 'omitnan'), median(Std(:), 'omitnan'));

%% 4. Source detection -----------------------------------------------------
Log('[4/8] Detecting sources (%g sigma, smoothing %g pix)', Args.DetectNsigma, Args.SmoothSigma);
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
Log('      %d sources detected', Nsrc);

%% 5. Flag saturated and blended stars ------------------------------------
IsSat   = Peak >= Args.NonLinFrac * SatLevel;
IsBlend = (Npeaks > 1) | (NNdist < Args.MinSep);
IsGood  = ~IsSat & ~IsBlend;
Log('[5/8] Flagging stars: good %d, saturated %d, blended %d', ...
    nnz(IsGood), nnz(IsSat), nnz(IsBlend & ~IsSat));

%% 6. Mask bad stars and fill with background + noise ---------------------
Log('[6/8] Masking saturated and blended stars');
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

% image used by moments/bank: same as ImForACF, NaNs replaced by background
BackLevel = median(Back(:), 'omitnan');
ImFilled  = ImForACF;
Bad = ~isfinite(ImFilled);
if any(Bad(:))
    if Args.FlattenBack
        ImFilled(Bad) = BackLevel;
    else
        ImFilled(Bad) = Back(Bad);
    end
end
Log('      %d pixels masked (%.2f%% of the image)', Nmask, 100*Nmask/numel(Mask));

%% 7. Regions and FWHM measurement -----------------------------------------
[Sec, GridSize] = local_make_sections(Nx, Ny, Args.CCDSEC, Args.Grid);
Nsec = size(Sec, 1);
Log('[7/8] Measuring FWHM in %d region(s) with %d method(s)', Nsec, Nm);

[Xc, Yc, MaskedFrac]  = deal(nan(Nsec, 1));
[Ngood, Nsat, Nblend] = deal(zeros(Nsec, 1));
LowStars              = false(Nsec, 1);
FWHMm   = nan(Nsec, Nm);
StatusM = false(Nsec, Nm);
NusedM  = nan(Nsec, Nm);

for k = 1:Nsec
    S = Sec(k,:);
    Xc(k) = (S(1) + S(2)) / 2;
    Yc(k) = (S(3) + S(4)) / 2;
    Ry = S(3):S(4);
    Rx = S(1):S(2);

    InSec = XY(:,1) >= S(1) & XY(:,1) <= S(2) & XY(:,2) >= S(3) & XY(:,2) <= S(4);
    Ngood(k)  = nnz(InSec & IsGood);
    Nsat(k)   = nnz(InSec & IsSat);
    Nblend(k) = nnz(InSec & IsBlend & ~IsSat);
    MaskedFrac(k) = mean(Mask(Ry, Rx), 'all');
    LowStars(k)   = Ngood(k) < Args.MinStars;

    for m = 1:Nm
        try
            switch Methods(m)
                case "acf"
                    RegSize = [numel(Ry), numel(Rx)];
                    MaxRad  = min(Args.MaxRadius, floor(min(RegSize)/2) - 1);
                    if MaxRad < 3*Args.Step
                        warning('measure_wholeImageFWHM_from_fits:SmallRegion', ...
                            'Region %d is too small for the ACF method; skipped.', k);
                        Log('      region %d/%d %-7s skipped (region too small)', k, Nsec, Methods(m));
                        continue
                    end
                    [F, Ns, Info] = AcfFun(ImForACF, 'CCDSEC', S, ...
                        'CorrFrac', Args.CorrFrac, 'Nsigma0', Args.Nsigma0, ...
                        'MaxRadius', MaxRad, 'Step', Args.Step, ...
                        'SatLevel', Args.MedianSatLevel);

                case "moments"
                    [F, Ns, Info] = imUtil.psf.fwhm_fromMoments(ImFilled, 'CCDSEC', S, ...
                        'MinSN', Args.MomMinSN, 'MomRadius', Args.MomRadius);

                case "bank"
                    BankArgs = {'CCDSEC', S, 'MinSN', Args.BankMinSN, ...
                        'SigmaVec', Args.BankSigmaVec, 'MinStars', Args.BankMinStars, ...
                        'PixScale', 1, 'Method', char(Args.BankMethod), ...
                        'MaxIter', Args.BankMaxIter};
                    if Args.BankUseBackMap
                        if Args.FlattenBack
                            BackArg = single(BackLevel);
                        else
                            BackArg = single(Back(Ry, Rx));
                        end
                        BankArgs = [BankArgs, {'Background', BackArg, ...
                            'Variance', single(Std(Ry, Rx).^2)}]; %#ok<AGROW>
                    end
                    [F, Ns, Info] = imUtil.psf.fwhm_fromBank(ImFilled, BankArgs{:});
            end

            FWHMm(k,m)  = F;
            NusedM(k,m) = Ns;
            if isstruct(Info) && isfield(Info, 'Status')
                St = logical(Info.Status);
            else
                St = true;               % fwhm_fromMoments has no Status field
            end
            StatusM(k,m) = St && isfinite(F);
            Log('      region %d/%d [X %d-%d, Y %d-%d] %-7s FWHM = %6.2f pix  (%s)', ...
                k, Nsec, S(1), S(2), S(3), S(4), Methods(m), F, ...
                local_tf(StatusM(k,m), 'ok', 'failed'));
        catch ME
            warning('measure_wholeImageFWHM_from_fits:MethodFailed', ...
                'Region %d, method %s: %s', k, Methods(m), ME.message);
            Log('      region %d/%d %-7s error: %s', k, Nsec, Methods(m), ME.message);
        end
    end
end

Result = table((1:Nsec).', Sec(:,1), Sec(:,2), Sec(:,3), Sec(:,4), Xc, Yc, ...
    FWHMm(:,1), StatusM(:,1), Ngood, Nsat, Nblend, MaskedFrac, LowStars, ...
    'VariableNames', {'Region','Xmin','Xmax','Ymin','Ymax','Xcenter','Ycenter', ...
    'FWHM','Status','NgoodStars','NsatStars','NblendStars','MaskedFrac','LowStars'});
for m = 1:Nm
    Result.("FWHM_"   + Methods(m)) = FWHMm(:,m);
    Result.("Status_" + Methods(m)) = StatusM(:,m);
    Result.("Nused_"  + Methods(m)) = NusedM(:,m);
end
if ~isempty(Args.PixScale)
    Result.FWHM_arcsec = Result.FWHM .* Args.PixScale;
    for m = 1:Nm
        Result.("FWHM_" + Methods(m) + "_arcsec") = FWHMm(:,m) .* Args.PixScale;
    end
end

Sources = table(XY(:,1), XY(:,2), Area, Peak, Npeaks, NNdist, IsSat, IsBlend, IsGood, ...
    'VariableNames', {'X','Y','Area','Peak','Npeaks','NNdist','IsSat','IsBlend','IsGood'});

Aux = struct('Sources', Sources, 'Mask', Mask, 'Back', Back, 'Std', Std, ...
    'CleanImage', Clean, 'SatLevel', SatLevel, 'Methods', Methods, ...
    'ACFEngine', ACFEngine, 'Sections', Sec);

%% 8. Report and plots ------------------------------------------------------
Log('[8/8] Reporting results');
if Args.Verbose
    fprintf('\n%s  [%d x %d pix]\n', Name, Nx, Ny);
    fprintf('Methods: %s (primary: %s) | ACF engine: %s | SatLevel: %g\n', ...
        strjoin(Methods, ', '), Methods(1), ACFEngine, SatLevel);
    fprintf('Sources: %d (good %d, saturated %d, blended %d) | masked pixels: %.2f%%\n\n', ...
        Nsrc, nnz(IsGood), nnz(IsSat), nnz(IsBlend & ~IsSat), 100*Nmask/numel(Mask));
    disp(Result)
end

if Args.Plot
    Log('      plotting figures');
    local_plot(Image, Sec, Result, Sources, GridSize, Name, Methods, Args.MapFontSize);
end
Log('measure_wholeImageFWHM_from_fits: done');
end


function local_log(Verbose, T0, Fmt, varargin)
% Print one stage message with the elapsed time since the start of the run.
if Verbose
    fprintf('[%7.2f s] %s\n', toc(T0), sprintf(Fmt, varargin{:}));
end
end


function Out = local_tf(Cond, A, B)
if Cond, Out = A; else, Out = B; end
end


%% ========================================================================
%  Local functions
%% ========================================================================

function [Image, SatHdr] = local_read_fits(FileName, HDU)
% Read a 2-D image HDU and look for a saturation keyword.
import matlab.io.*
if ~isfile(FileName)
    error('measure_wholeImageFWHM_from_fits:FileNotFound', 'File not found: %s', FileName);
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
        error('measure_wholeImageFWHM_from_fits:NoImage', 'No 2-D image HDU found in %s', FileName);
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
    error('measure_wholeImageFWHM_from_fits:RegionArgs', 'Specify either CCDSEC or Grid, not both.');
end

if ~isempty(CCDSEC)
    if size(CCDSEC, 2) ~= 4
        error('measure_wholeImageFWHM_from_fits:CCDSEC', 'CCDSEC must be N-by-4: [Xmin Xmax Ymin Ymax].');
    end
    Sec = round(double(CCDSEC));
    Orig = Sec;
    Sec(:,[1 2]) = min(max(Sec(:,[1 2]), 1), Nx);
    Sec(:,[3 4]) = min(max(Sec(:,[3 4]), 1), Ny);
    if ~isequal(Sec, Orig)
        warning('measure_wholeImageFWHM_from_fits:CCDSECclipped', 'Some CCDSEC rows were clipped to the image bounds.');
    end
    if any(Sec(:,2) <= Sec(:,1) | Sec(:,4) <= Sec(:,3))
        error('measure_wholeImageFWHM_from_fits:CCDSEC', 'Each CCDSEC row needs Xmin<Xmax and Ymin<Ymax.');
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



function local_plot(Image, Sec, Result, Sources, GridSize, Name, Methods, MapFontSize)
% Image with sources and regions, FWHM map for multi-region runs, and a
% method comparison when several methods were run.
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
title(sprintf('%s  (FWHM per region, method: %s)', Name, Methods(1)), 'Interpreter', 'none');
hold off

if height(Result) > 1
    % one FWHM map per chosen method, side by side in one figure
    Nm  = numel(Methods);
    Nc  = min(Nm, 3);
    Nr  = ceil(Nm / Nc);
    figure('Name', 'FWHM maps', 'Color', 'w');
    TL = tiledlayout(Nr, Nc, 'TileSpacing', 'compact', 'Padding', 'compact');
    for m = 1:Nm
        ax = nexttile(TL);
        Fm = Result.("FWHM_" + Methods(m));
        if ~isempty(GridSize)
            Map = reshape(Fm, GridSize(1), GridSize(2)).';
            Xv  = Result.Xcenter(1:GridSize(1));
            Yv  = Result.Ycenter(1:GridSize(1):end);
            imagesc(ax, Xv, Yv, Map, 'AlphaData', isfinite(Map));
            axis(ax, 'image'); set(ax, 'YDir', 'normal', 'Color', [0.6 0.75 0.9]);
            colormap(ax, hot); colorbar(ax);
            CL = clim(ax);
            for k = 1:height(Result)
                if Methods(m) == "acf"
                    Lbl = sprintf('%.2f\n(%d)', Fm(k), Result.NgoodStars(k));
                else
                    Lbl = sprintf('%.2f\n(%d)', Fm(k), Result.("Nused_" + Methods(m))(k));
                end
                % 'hot' is bright at the top: black text on light cells, white on dark
                if isfinite(Fm(k)) && diff(CL) > 0 && (Fm(k) - CL(1))/diff(CL) > 0.55
                    TxtCol = 'k';
                else
                    TxtCol = 'w';
                end
                text(ax, Result.Xcenter(k), Result.Ycenter(k), Lbl, ...
                    'HorizontalAlignment', 'center', 'Color', TxtCol, ...
                    'FontSize', MapFontSize, 'FontWeight', 'bold');
            end
        else
            scatter(ax, Result.Xcenter, Result.Ycenter, 200, Fm, 'filled', 'MarkerEdgeColor', 'k');
            axis(ax, 'equal'); colormap(ax, hot); colorbar(ax); grid(ax, 'on');
            for k = 1:height(Result)
                text(ax, Result.Xcenter(k), Result.Ycenter(k), sprintf('  %.2f', Fm(k)), ...
                    'FontSize', MapFontSize);
            end
        end
        xlabel(ax, 'X [pix]'); ylabel(ax, 'Y [pix]');
        if Methods(m) == "acf"
            title(ax, 'FWHM map: acf (good stars)');
        else
            title(ax, sprintf('FWHM map: %s (stars used)', Methods(m)));
        end
    end
    title(TL, sprintf('%s', Name), 'Interpreter', 'none');
end

if numel(Methods) > 1
    figure('Name', 'FWHM: method comparison', 'Color', 'w');
    hold on
    for m = 1:numel(Methods)
        plot(Result.Region, Result.("FWHM_" + Methods(m)), '-o', ...
            'LineWidth', 1.2, 'DisplayName', Methods(m));
    end
    hold off
    grid on; legend('Location', 'best');
    xlabel('Region'); ylabel('FWHM [pix]');
    title('FWHM per region for each method');
    if height(Result) > 1
        xticks(Result.Region);
    end
end
end