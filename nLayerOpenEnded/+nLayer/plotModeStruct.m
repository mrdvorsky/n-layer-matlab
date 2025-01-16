function [] = plotModeStruct(modeStruct, options)
%PLOTMODESTRUCT Plot the fields defined by the input modeStruct.
%   Detailed explanation goes here

arguments
    modeStruct(1, 1);

    options.SizeX(1, 1) {mustBePositive};
    options.SizeY(1, 1) {mustBePositive};
    options.NumPointsX(1, 1) {mustBePositive, mustBeInteger} = 999;
    options.NumPointsY(1, 1) {mustBePositive, mustBeInteger} = 999;
end

%% Check Size
if ~isfield(options, "SizeX")
    options.SizeX = 1.1 * max([modeStruct.ApertureWidth]);
end
if ~isfield(options, "SizeY")
    options.SizeY = 1.1 * max([modeStruct.ApertureWidth]);
end

%% Calculate x, y, kx, and ky
x(:, 1) = 1.2 * options.SizeX * linspace(-0.5, 0.5, options.NumPointsX);
y(1, :) = 1.2 * options.SizeY * linspace(-0.5, 0.5, options.NumPointsY);

[kx, ky] = fftCoordinates(x, y, ApplyFftShift=true);
kr = hypot(kx, ky);
kphi = atan2(ky, kx);

%% Extract Modes
modeLabel = modeStruct.ModeLabel;
ExHat = modeStruct.ExSpec(kx, ky, kr, kphi);
EyHat = modeStruct.EySpec(kx, ky, kr, kphi);
HertzHat = modeStruct.EzSpec(kx, ky, kr, kphi);

scaleFactor = numel(ExHat) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
    ./ (2*pi);
Ex = fftshift(ifft2(ifftshift(ExHat))) * scaleFactor;
Ey = fftshift(ifft2(ifftshift(EyHat))) * scaleFactor;
Ez = fftshift(ifft2(ifftshift(HertzHat))) * scaleFactor;

% Check symmetry.
symX1 = "None";
if rms(abs(Ex + flip(Ex, 2)), "all") < 0.05
    symX1 = "PEC";
elseif rms(abs(Ex - flip(Ex, 2)), "all") < 0.05
    symX1 = "PMC";
end

symY1 = "None";
if rms(abs(Ey + flip(Ey, 1)), "all") < 0.05
    symY1 = "PEC";
elseif rms(abs(Ey - flip(Ey, 1)), "all") < 0.05
    symY1 = "PMC";
end

modeAmp = (sum(ExHat.*flip(flip(ExHat, 1), 2), "all") ...
    + sum(EyHat.*flip(flip(EyHat, 1), 2), "all")) ...
    .* (abs(kx(2) - kx(1)) .* abs(ky(2) - ky(1)));

% modeAmp = sum(Ex(:).^2 + Ey(:).^2) ...
%     .* (abs(x(2) - x(1)) .* abs(y(2) - y(1)));


% subplot(1, 2, 1);
plotVectorField(x, y, Ex, Ey);

% subplot(1, 2, 2);
% showImage(x, y, Ez, DisplayFormat="Imag");

title(sprintf("%s: (%.3f) (%s, %s)", modeLabel, modeAmp, symX1, symY1));


end

