function [] = plotModeStruct(modeStruct, options)
%PLOTMODESTRUCT Plot the fields defined by the input modeStruct.
%   Detailed explanation goes here

arguments
    modeStruct(1, 1);

    options.SizeX(1, 1) {mustBeReal} = 1;
    options.SizeY(1, 1) {mustBeReal} = 1;
    options.NumPointsX(1, 1) {mustBePositive, mustBeInteger} = 1000;
    options.NumPointsY(1, 1) {mustBePositive, mustBeInteger} = 1000;

    options.PlotCoordinates {mustBeMember(options.PlotCoordinates, ...
        ["Rectangular", "Polar"])} = "Rectangular";
end

%% Calculate x, y, kx, and ky
x(:, 1) = 1.2 * options.SizeX * linspace(-0.5, 0.5, options.NumPointsX);
y(1, :) = 1.2 * options.SizeY * linspace(-0.5, 0.5, options.NumPointsY);

[kx, ky] = fftCoordinates(x, y);
kRho = hypot(kx, ky);
kPhi = atan2(ky, kx);

%% Extract Modes
specEx_All = [modeStruct.SpecEx_TE; modeStruct.SpecEx_TM; modeStruct.SpecEx_Hybrid];
specEy_All = [modeStruct.SpecEy_TE; modeStruct.SpecEy_TM; modeStruct.SpecEy_Hybrid];
modeScale_All = [modeStruct.PhaseScaleFactor_TE; modeStruct.PhaseScaleFactor_TM; modeStruct.PhaseScaleFactor_Hybrid];
modeLabels_All = [modeStruct.ModeLabels_TE; modeStruct.ModeLabels_TM; modeStruct.ModeLabels_Hybrid];
coordLabels = ["x", "y"];
for ii = 1:numel(specEx_All)
    ExHat = specEx_All{ii}(kx, ky, kRho, kPhi) + zeros(size(kRho));
    EyHat = specEy_All{ii}(kx, ky, kRho, kPhi) + zeros(size(kRho));

    scaleFactor = numel(ExHat) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
        ./ (2*pi) .* modeScale_All(ii);
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

    modeAmp = trapz(x, trapz(y, Ex.^2 + Ey.^2, 2), 1);

    % Plot Fields
    figure;
    subplot(2, 2, 1);
    showImage(x, y, Ex, DisplayFormat="Real");
    colormap colormapPlusMinus;
    clim(colorscale * [-1, 1]);
    title(sprintf("E_{%s} for %s (%.3f)", coordLabels(1), modeLabels_All(ii), modeAmp));

    subplot(2, 2, 3);
    showImage(x, y, Ex, DisplayFormat="Imag");
    colormap colormapPlusMinus;
    clim(colorscale * [-1, 1]);
    title(sprintf("E_{%s} for %s (%.3f)", coordLabels(1), modeLabels_All(ii), modeAmp));

    subplot(2, 2, 2);
    showImage(x, y, Ey, DisplayFormat="Real");
    colormap colormapPlusMinus;
    clim(colorscale * [-1, 1]);
    title(sprintf("E_{%s} for %s (%.3f)", coordLabels(2), modeLabels_All(ii), modeAmp));

    subplot(2, 2, 4);
    showImage(x, y, Ey, DisplayFormat="Imag");
    colormap colormapPlusMinus;
    clim(colorscale * [-1, 1]);
    title(sprintf("E_{%s} for %s (%.3f)", coordLabels(2), modeLabels_All(ii), modeAmp));
end




end

