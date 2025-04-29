clc;
clear;
close all;

%% Inputs
m = 2;
n = 2;

NL = nLayerRectangular(m, n, waveguideBand="Ka", modeSymmetryX="None", modeSymmetryY="None");
wgA = NL.waveguideA;
wgB = NL.waveguideB;

mode = NL.modeStructs(strcmp(NL.modeLabels, "TE_{1,1}"));

kc = mode.CutoffWavenumber;

kx(:, 1) = 160 * linspace(-1, 1, 10001);
ky(1, :) = 160 * linspace(-1, 1, 10001);

kr = hypot(kx, ky);
kphi = angle(kx + 1j*ky);

%% Analytical
spec1 = mode.ExSpec(kx, ky) + 0 .* mode.EySpec(kx, ky);

% spec1 = 

%% Spectrums Area
% spec2_kx = reshape(...
%     nufftn(mode.Hertz_Hz .* mode.Hertz_w .* cos(mode.Hertz_ang), [mode.Hertz_x, mode.Hertz_y], ...
%     {kx(:) ./ (2*pi), ky(:) ./ (2*pi)}), ...
%     numel(kx), numel(ky)) ...
%     .* cos(kphi);
% 
% spec2_ky = reshape(...
%     nufftn(mode.Hertz_Hz .* mode.Hertz_w .* sin(mode.Hertz_ang), [mode.Hertz_x, mode.Hertz_y], ...
%     {kx(:) ./ (2*pi), ky(:) ./ (2*pi)}), ...
%     numel(kx), numel(ky)) ...
%     .* sin(kphi);
% 
% spec2 = -1j * (spec2_kx + spec2_ky);

%% Spatial
[x, y] = fftCoordinates(kx, ky, ApplyFftShift=true);

scale = abs(kx(2) - kx(1)) * abs(ky(2) - ky(1)) .* numel(spec1);
spat1 = fftshift(ifft2(ifftshift(spec1))) .* scale;

%% Plot Spectrum
figure;
showImage(kx, ky, spec1, DisplayFormat="Magnitude");
colormap jet;
grid on;

% figure;
% showImage(kx, ky, spec2, DisplayFormat="Magnitude");
% colormap jet;
% grid on;
% 
% figure;
% showImage(kx, ky, (spec1 - spec2) ./ max(abs(spec1(:))), DisplayFormat="dB");
% colormap jet;
% grid on;
% clim([-300, 0]);

%% Plot Spatial
xPlotVals = (x < wgA) & (x > -wgA);
yPlotVals = (y < wgB) & (y > -wgB);

figure;
showImage(x(xPlotVals), y(yPlotVals), spat1(xPlotVals, yPlotVals), DisplayFormat="Magnitude");
colormap jet;
grid on;




