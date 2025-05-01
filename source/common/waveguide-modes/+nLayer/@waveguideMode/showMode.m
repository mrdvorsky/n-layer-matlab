function [] = showMode(O, options)
%Plot the fields defined by this "waveguideMode" object.
% This function simply plots the tangential fields over the waveguide
% aperture.
%
% Author: Matt Dvorsky

arguments
    O nLayer.waveguideMode;

    options.SizeX(1, 1) {mustBePositive};
    options.SizeY(1, 1) {mustBePositive};
    options.NumPointsX(1, 1) {mustBePositive, mustBeInteger} = 500;
    options.NumPointsY(1, 1) {mustBePositive, mustBeInteger} = 500;

    options.ArrowDecimationFactorX(1, 1) {mustBePositive, mustBeInteger} = 10;
    options.ArrowDecimationFactorY(1, 1) {mustBePositive, mustBeInteger} = 10;
end

%% Check Size
if ~isfield(options, "SizeX")
    options.SizeX = 1.1 * O.apertureSize;
end
if ~isfield(options, "SizeY")
    options.SizeY = 1.1 * O.apertureSize;
end

%% Calculate x, y, kx, and ky
x(:, 1) = options.SizeX * linspace(-0.5, 0.5, options.NumPointsX);
y(1, :) = options.SizeY * linspace(-0.5, 0.5, options.NumPointsY);

[kx, ky] = fftCoordinates(x, y, ApplyFftShift=true);
kr = hypot(kx, ky);
kphi = atan2(ky, kx);

%% Calculate Fields
ExHat = O.ExSpec(kx, ky, kr, kphi);
EyHat = O.EySpec(kx, ky, kr, kphi);

intScaleFactor = numel(kr) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
    ./ (2*pi);

Ex = fftshift(ifft2(ifftshift(ExHat))) * intScaleFactor;
Ey = fftshift(ifft2(ifftshift(EyHat))) * intScaleFactor;

%% Plot Field Magnitudes
Er = hypot(Ex, Ey);
showImage(x, y, Er, DisplayFormat="Magnitude");

%% Plot Vectors
dx = abs(x(2) - x(1));
dy = abs(y(2) - y(1));
decX = options.ArrowDecimationFactorX;
decY = options.ArrowDecimationFactorY;

arrowScaleFactor = 0.7*max(dx*decX, dy*decY) ./ max(Er(:));
uq = arrowScaleFactor * real(Ex);
vq = arrowScaleFactor * real(Ey);

xq = x - 0.5*uq;
yq = y - 0.5*vq;

xq = xq(1:decX:end, 1:decY:end);
yq = yq(1:decX:end, 1:decY:end);
vq = vq(1:decX:end, 1:decY:end);
uq = uq(1:decX:end, 1:decY:end);

hold on;
quiver(xq(:), yq(:), uq(:), vq(:), "off", "k");

end

