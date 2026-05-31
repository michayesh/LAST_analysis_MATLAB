function plotPolarIntensityMap(radii, azimuths_deg, intensity)
% PLOTPOLARINTENSITYMAP Plots a 2D intensity map in polar coordinates.
%
% Inputs:
%   radii         - Vector of radial distances (length M)
%   azimuths_deg  - Vector of azimuth angles in degrees (length N)
%   intensity     - 2D matrix of intensity values (M x N),
%                   rows correspond to radii, columns to azimuths
%
% Example:
%   r = linspace(0, 1, 50);
%   theta = linspace(0, 360, 72);
%   [R, T] = meshgrid(r, theta);
%   I = R' .* cosd(T').^2;
%   plotPolarIntensityMap(r, theta, I);

    %% --- Input validation ---
    if ~isvector(radii) || ~isvector(azimuths_deg)
        error('radii and azimuths_deg must be vectors.');
    end
    radii        = radii(:);           % ensure column vector
    azimuths_deg = azimuths_deg(:)';   % ensure row vector

    [M, N] = size(intensity);
    if M ~= numel(radii) || N ~= numel(azimuths_deg)
        error(['intensity must be of size [numel(radii) x numel(azimuths_deg)] ' ...
               '= [%d x %d], but got [%d x %d].'], ...
               numel(radii), numel(azimuths_deg), M, N);
    end

    %% --- Convert polar grid to Cartesian ---
    azimuths_rad = deg2rad(azimuths_deg);          % 1 x N

    % Build 2-D grids  (M x N each)
    [AZ, R] = meshgrid(azimuths_rad, radii);

    X = R .* cos(AZ);
    Y = R .* sin(AZ);

    %% --- Plot ---
    figure('Name', 'Polar Intensity Map', 'Color', 'k');

    % pcolor needs one extra row/column to define patch edges;
    % use surf(..., 'EdgeColor','none') on the interior grid instead.
    surf(X, Y, intensity, 'EdgeColor', 'none');
    view(2);                   % top-down view → 2-D appearance

    axis equal tight;
    colormap(gca, 'jet');
    cb = colorbar;
    cb.Color = 'w';
    cb.Label.String = 'Intensity';
    cb.Label.Color  = 'w';

    xlabel('X', 'Color', 'w');
    ylabel('Y', 'Color', 'w');
    title('2D Polar Intensity Map', 'Color', 'w');

    ax = gca;
    ax.Color        = 'k';
    ax.XColor       = 'w';
    ax.YColor       = 'w';
    ax.GridColor    = 'w';
    ax.GridAlpha    = 0.2;
    grid on;
end