clc;
clear;
close all;

%% Inputs
krMax = 1000;

phiR = linspace(0, pi/2, 1001);
phiI = linspace(-1, 1, 101);

%% nLayer
wgA = 1;
wgB = 1/2;

%% Get Spectrum Functions
[spEx1, spEy1] = nLayerOpenEnded.getSpectrumRectangular(...
    wgA, wgB, 3, 2, "TM");
[spEx2, spEy2] = nLayerOpenEnded.getSpectrumRectangular(...
    wgA, wgB, 3, 0, "TM");

kr(:, 1) = krMax;
kx = kr .* cos(phiR);
ky = kr .* sin(phiR);

spec = spEy1(kx, ky);

%% Plot Spectrum
figure;
plot(rad2deg(phiR), real(spec), "", LineWidth=1.5)
grid on;







