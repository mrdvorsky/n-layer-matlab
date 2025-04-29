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
    options.NumPointsX(1, 1) {mustBePositive, mustBeInteger} = 625;
    options.NumPointsY(1, 1) {mustBePositive, mustBeInteger} = 625;
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

scaleFactor = numel(ExHat) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
    ./ (2*pi);
Ex = fftshift(ifft2(ifftshift(ExHat))) * scaleFactor;
Ey = fftshift(ifft2(ifftshift(EyHat))) * scaleFactor;

%% Plot Fields
showImage(x, y, sqrt(Ex.^2 + Ey.^2), DisplayFormat="Magnitude");

end

