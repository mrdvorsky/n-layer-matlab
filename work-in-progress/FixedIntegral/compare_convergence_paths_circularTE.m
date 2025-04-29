clc;
clear;
close all;

%% Inputs
Nc = 11;
L = 3;
Lc = 2;

Lch = Lc;
Lcw = 20;


Nrho = 16000;
Nphi = 50;

%% Structure
er = 4 - 0.004j;
ur = 1;
thk = 1.2;
f = 40;
k0 = f * (2*pi/299.792458);

%% Contour Paths
krToKrc = @(x) x + 1j * (Lch.*Lcw) * x ./ sqrt((Lch + x.^2).*(3*Lcw.^2 + x.^2));
krToKrcPrime = @(x) 1 + 1j * (Lch.*Lcw) * (3*Lch*Lcw.^2 - x.^4) ./ ((Lch + x.^2).*(3*Lcw.^2 + x.^2)).^1.5;

%% nLayer
wgR = 30/128 * 25.4;

%% Get Spectrum Functions
[specX1, specY1] = nLayer.getSpectrumCircular(...
    wgR, 0, 1, "TE");

%% Integral Old
tic;
Ih = 2*pi * integral(@(x) gamH(x, k0, er, ur, thk) .* ...
    x .* sum(specY1(x, 0, x, 0).^2, 3), ...
    0, inf, RelTol=1e-7);
toc;

%% Integral 1
[krc1, krc_weights_H1] = getContourWeights_circularTE(...
    Nc, Nrho, Nphi, L, Lc, wgR, krToKrc, krToKrcPrime);

I1h = sum(krc_weights_H1 .* gamH(krc1, k0, er, ur, thk));

%% Error
errH1 = db(I1h - Ih)

%% Plot Gamma
[krPlot] = fejer2_halfOpen(1000, L);
krcPlot = krPlot + Lc .* 1j .* krPlot ./ (Lc + krPlot);

figure;
plot(1:numel(krcPlot), real(gamH(krcPlot, k0, er, ur, thk) ./ sqrt(1 + krcPlot.^2)), "", LineWidth=1.5);
hold on;
plot(1:numel(krcPlot), imag(gamH(krcPlot, k0, er, ur, thk) ./ sqrt(1 + krcPlot.^2)), "", LineWidth=1.5);
grid on;
title("gamH");

figure;
plot(1:numel(krcPlot), real(gamE(krcPlot, k0, er, ur, thk) .* sqrt(1 + krcPlot.^2)), "", LineWidth=1.5);
hold on;
plot(1:numel(krcPlot), imag(gamE(krcPlot, k0, er, ur, thk) .* sqrt(1 + krcPlot.^2)), "", LineWidth=1.5);
grid on;
title("gamE");

%% Plotting
figure;
plot(1:numel(krc1), real(krc_weights_H1), "", LineWidth=1.5);
hold on;
plot(1:numel(krc1), imag(krc_weights_H1), "", LineWidth=1.5);grid on;
title("gamH Weights");





