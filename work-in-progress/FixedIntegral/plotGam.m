clc;
clear;
close all;

%% Inputs
Nc = 50;
L = 1;
Lc = 1.0;

Nrho = 16000;
Nphi = 450;

%% Structure
er = 2 - 0.0001j;
ur = 1;
thk = 2.5;
f = 30;
k0 = f * (2*pi/299.792458);

%% nLayer
wgA = 1;
wgB = 0.5;

%% Get Spectrum Functions
[specX1, specY1] = nLayerOpenEnded.getSpectrumRectangular(...
    wgA, wgB, 1, 0, "TE");

%% Integral Accurate
[krc, krc_weights_E1, krc_weights_H1] = getContourWeights(...
    Nc, Nrho, Nphi, L, Lc, wgA, wgB);

%% Plotting
figure;
plot(1:numel(krc), real(gamH(krc, k0, er, ur, thk) ./ sqrt(1 + krc.^2)), "", LineWidth=1.5);
hold on;
plot(1:numel(krc), imag(gamH(krc, k0, er, ur, thk) ./ sqrt(1 + krc.^2)), "", LineWidth=1.5);
grid on;
title("gamH");

figure;
plot(1:numel(krc), real(gamE(krc, k0, er, ur, thk) .* sqrt(1 + krc.^2)), "", LineWidth=1.5);
hold on;
plot(1:numel(krc), imag(gamE(krc, k0, er, ur, thk) .* sqrt(1 + krc.^2)), "", LineWidth=1.5);
grid on;
title("gamE");







