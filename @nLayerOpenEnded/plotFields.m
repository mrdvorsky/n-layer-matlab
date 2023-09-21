function [] = plotFields(modeStruct, x, y)
%PLOTFIELDS Plot the fields defined by modeStruct.
%   Detailed explanation goes here

arguments
    modeStruct(1, 1);
    x(:, 1) {mustBeReal};
    y(1, :) {mustBeReal};
end

%% Calculate kx and ky
[kx, ky] = fftCoordinates(x, y);
kRho = hypot(kx, ky);
kPhi = atan2(ky, kx);

%% Extract Modes
specEx_All = [modeStruct.SpecEx_TE; modeStruct.SpecEx_TM; modeStruct.SpecEx_Hybrid];
specEy_All = [modeStruct.SpecEy_TE; modeStruct.SpecEy_TM; modeStruct.SpecEy_Hybrid];
modeLabels_All = [modeStruct.ModeLabels_TE; modeStruct.ModeLabels_TM; modeStruct.ModeLabels_Hybrid];
for ii = 1:numel(specEx_All)
    ExHat = specEx_All{ii}(kx, ky, kRho, kPhi) + zeros(size(kRho));
    EyHat = specEy_All{ii}(kx, ky, kRho, kPhi) + zeros(size(kRho));

    scaleFactor = numel(ExHat) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
        ./ (2*pi);
    Ex = ifftshift(ifft2(ExHat)) * scaleFactor;
    Ey = ifftshift(ifft2(EyHat)) * scaleFactor;

    modeAmp = trapz(x, trapz(y, abs(Ex).^2 + abs(Ey).^2, 2), 1);

    figure;
    subplot(2, 2, 1);
    showImage(x, y, Ex, DisplayFormat="Real");
    colormap jet;
    title(sprintf("E_x for %s (%.3f)", modeLabels_All(ii), modeAmp));

    subplot(2, 2, 3);
    showImage(x, y, Ex, DisplayFormat="Imag");
    colormap jet;
    title(sprintf("E_x for %s (%.3f)", modeLabels_All(ii), modeAmp));

    subplot(2, 2, 2);
    showImage(x, y, Ey, DisplayFormat="Real");
    colormap jet;
    title(sprintf("E_y for %s (%.3f)", modeLabels_All(ii), modeAmp));

    subplot(2, 2, 4);
    showImage(x, y, Ey, DisplayFormat="Imag");
    colormap jet;
    title(sprintf("E_y for %s (%.3f)", modeLabels_All(ii), modeAmp));
end




end

