clc;
clear;
close all;

%% Inputs
m = 2;
n = 2;

NL = nLayerRectangular(m, n, waveguideBand="Ka", modeSymmetryX="None", modeSymmetryY="None");
wgA = NL.waveguideA;
wgB = NL.waveguideB;

mode = NL.modeStructs(end);     % Last mode is TMmn.

kc = mode.CutoffWavenumber;

kx(:, 1) = 1080 * linspace(-1, 1, 1000);
ky(1, :) = 1080 * linspace(-1, 1, 1000);

kr = hypot(kx, ky);
kphi = angle(kx + 1j*ky);

%% Analytical
spec1 = cos(kphi) .* mode.ExSpec(kx, ky) + sin(kphi) .* mode.EySpec(kx, ky);

%% Spectrums Area
spec2 = reshape(...
    nufftn(mode.Hertz_Ez .* mode.Hertz_w, [mode.Hertz_x, mode.Hertz_y], ...
    {kx(:) ./ (2*pi), ky(:) ./ (2*pi)}), ...
    numel(kx), numel(ky)) ...
    .* kr ./ (kr.^2 - kc.^2);

%% Plot Spectrum
figure;
showImage(kx, ky, spec1, DisplayFormat="Magnitude");
colormap jet;
grid on;

figure;
showImage(kx, ky, spec2, DisplayFormat="Magnitude");
colormap jet;
grid on;

figure;
showImage(kx, ky, (spec1 - spec2) ./ max(abs(spec1(:))), DisplayFormat="dB");
colormap jet;
grid on;
clim([-300, 0]);





