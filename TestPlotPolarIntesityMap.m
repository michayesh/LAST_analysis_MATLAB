r     = linspace(0, 1, 60);
theta = linspace(0, 360, 120);

[R, T] = meshgrid(r, theta);
I = R' .* (1 + cosd(2*T'))./ 2;   % donut-like pattern

plotPolarIntensityMap(r, theta, I);