function [] = showMode(self, options)
%Plot the fields defined by this "waveguideMode" object.
% This function simply plots the tangential fields over the waveguide
% aperture.
%
% Author: Matt Dvorsky

arguments
    self waveguideMode;

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

%% Calculate x and y
x(:, 1) = options.SizeX * linspace(-0.5, 0.5, options.NumPointsX);
y(1, :) = options.SizeY * linspace(-0.5, 0.5, options.NumPointsY);

%% Calculate Fields
[Ex, Ey] = self.getModeFields(x, y);

%% Plot Field Magnitudes
Er = hypot(Ex, Ey);
showImage(x, y, Er, DisplayFormat="Magnitude", ...
    Axis=options.Axis, ...
    Interpolation="Bilinear");

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

xq = xq(round(0.5*decX):decX:end, round(0.5*decY):decY:end);
yq = yq(round(0.5*decX):decX:end, round(0.5*decY):decY:end);
vq = vq(round(0.5*decX):decX:end, round(0.5*decY):decY:end);
uq = uq(round(0.5*decX):decX:end, round(0.5*decY):decY:end);

hold(options.Axis, "on");
quiver(options.Axis, xq(:), yq(:), uq(:), vq(:), "off", "k");

end

