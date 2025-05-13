function [] = showMode(self, options)
%Plot the fields defined by this "waveguideMode" object.
% This function simply plots the tangential fields over the waveguide
% aperture.
%
% Author: Matt Dvorsky

arguments
    self nLayer.waveguideMode;

    options.SizeX(1, 1) {mustBePositive};
    options.SizeY(1, 1) {mustBePositive};
    options.NumPointsX(1, 1) {mustBePositive, mustBeInteger} = 500;
    options.NumPointsY(1, 1) {mustBePositive, mustBeInteger} = 500;

    options.ArrowDecimationFactorX(1, 1) {mustBePositive, mustBeInteger} = 10;
    options.ArrowDecimationFactorY(1, 1) {mustBePositive, mustBeInteger} = 10;

    options.Axis(1, 1) matlab.graphics.axis.Axes;
end

%% Check Inputs
if ~isfield(options, "Axis")
    options.Axis = gca();
end

if ~isfield(options, "SizeX")
    options.SizeX = 1.1 * self.apertureSize;
end
if ~isfield(options, "SizeY")
    options.SizeY = 1.1 * self.apertureSize;
end

%% Calculate x, y, kx, and ky
x(:, 1) = options.SizeX * linspace(-0.5, 0.5, options.NumPointsX);
y(1, :) = options.SizeY * linspace(-0.5, 0.5, options.NumPointsY);

[kx, ky] = fftCoordinates(x, y, ApplyFftShift=true);
kr = hypot(kx, ky);
kphi = atan2(ky, kx);

%% Calculate Fields
ExHat = self.ExSpec(kx, ky, kr, kphi);
EyHat = self.EySpec(kx, ky, kr, kphi);

intScaleFactor = numel(kr) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
    ./ (2*pi);

Ex = fftshift(ifft2(ifftshift(ExHat))) * intScaleFactor;
Ey = fftshift(ifft2(ifftshift(EyHat))) * intScaleFactor;

%% Plot Field Magnitudes
Er = hypot(Ex, Ey);
showImage(x, y, Er, DisplayFormat="Magnitude", ...
    Axis=options.Axis);

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

hold(options.Axis, "on");
quiver(options.Axis, xq(:), yq(:), uq(:), vq(:), "off", "k");

end

