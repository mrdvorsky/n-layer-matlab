function [] = plotModeStruct(modeStructs, options)
%PLOTMODESTRUCT Plot the fields defined by the input modeStruct.
%   Detailed explanation goes here

arguments
    modeStructs(:, 1);

    options.SizeX(1, 1) {mustBePositive};
    options.SizeY(1, 1) {mustBePositive};
    options.NumPointsX(1, 1) {mustBePositive, mustBeInteger} = 999;
    options.NumPointsY(1, 1) {mustBePositive, mustBeInteger} = 999;

    options.PlotCoordinates string {mustBeMember(options.PlotCoordinates, ...
        ["Rectangular", "Polar"])} = "Rectangular";
end

%% Check Size
if ~isfield(options, "SizeX")
    options.SizeX = 1.1 * max([modeStructs.ApertureWidth]);
end
if ~isfield(options, "SizeY")
    options.SizeY = 1.1 * max([modeStructs.ApertureWidth]);
end

%% Calculate x, y, kx, and ky
x(:, 1) = 1.2 * options.SizeX * linspace(-0.5, 0.5, options.NumPointsX);
y(1, :) = 1.2 * options.SizeY * linspace(-0.5, 0.5, options.NumPointsY);

[kx, ky] = fftCoordinates(x, y);
kr = hypot(kx, ky);
kphi = atan2(ky, kx);

%% Extract Modes
modeLabels = [modeStructs.ModeLabel];
coordLabels = ["x", "y"];
for ii = flip(1:numel(modeStructs))
    ExHat = modeStructs(ii).ExSpec(kx, ky, kr, kphi) + zeros(size(kr));
    EyHat = modeStructs(ii).EySpec(kx, ky, kr, kphi) + zeros(size(kr));

    scaleFactor = numel(ExHat) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
        ./ (2*pi);
    Ex = ifftshift(ifft2(ExHat)) * scaleFactor;
    Ey = ifftshift(ifft2(EyHat)) * scaleFactor;

    if strcmp(options.PlotCoordinates, "Polar")
        phi = angle(x + 1j*y);
        Er = Ex .* cos(phi) + Ey .* sin(phi);
        Ephi = -Ex .* sin(phi) + Ey .* cos(phi);

        Ex = Er;
        Ey = Ephi;

        coordLabels = ["\rho", "\phi"];
    end

    colorscale = max(abs(cat(3, real(Ex), imag(Ex), real(Ey), imag(Ey))), [], "all");

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

    symY2 = "None";
    if rms(abs(Ex + flip(Ex, 1)), "all") < 0.05
        symY2 = "PMC";
    elseif rms(abs(Ex - flip(Ex, 1)), "all") < 0.05
        symY2 = "PEC";
    end

    symX2 = "None";
    if rms(abs(Ey + flip(Ey, 2)), "all") < 0.05
        symX2 = "PMC";
    elseif rms(abs(Ey - flip(Ey, 2)), "all") < 0.05
        symX2 = "PEC";
    end

    % modeAmp = trapz(x, trapz(y, Ex.^2 + Ey.^2, 2), 1);
    % modeAmp = trapz(x, trapz(y, Ex.*rot90(Ex, 2) + Ey.*rot90(Ey, 2), 2), 1);
    modeAmp = sum(ExHat(:).^2 + EyHat(:).^2) ...
        .* (abs(kx(2) - kx(1)) .* abs(ky(2) - ky(1)));

    % % Plot Fields
    % figure;
    % subplot(2, 2, 1);
    % showImage(x, y, Ex, DisplayFormat="Real");
    % colormap colormapPlusMinus;
    % clim(colorscale * [-1, 1]);
    % title(sprintf("E_{%s} for %s (%.3f)", coordLabels(1), modeLabels(ii), modeAmp));
    % 
    % subplot(2, 2, 3);
    % showImage(x, y, Ex, DisplayFormat="Imag");
    % colormap colormapPlusMinus;
    % clim(colorscale * [-1, 1]);
    % title(sprintf("E_{%s} for %s (%.3f)", coordLabels(1), modeLabels(ii), modeAmp));
    % 
    % subplot(2, 2, 2);
    % showImage(x, y, Ey, DisplayFormat="Real");
    % colormap colormapPlusMinus;
    % clim(colorscale * [-1, 1]);
    % title(sprintf("E_{%s} for %s (%s, %s)", coordLabels(2), modeLabels(ii), symX1, symY1));
    % 
    % subplot(2, 2, 4);
    % showImage(x, y, Ey, DisplayFormat="Imag");
    % colormap colormapPlusMinus;
    % clim(colorscale * [-1, 1]);
    % title(sprintf("E_{%s} for %s (%s, %s)", coordLabels(2), modeLabels(ii), symX2, symY2));

    figure;
    plotVectorField(x, y, Ex, Ey);
    title(sprintf("%s: (%.3f) (%s, %s)", modeLabels(ii), modeAmp, symX1, symY1));

end




end

