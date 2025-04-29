clc;
clear;
close all;

%% Inputs
Nc = 50;
L = 2;
Lc = 1;

Nrho = 8000;
Nphi = 250;

%% Structure
er = 200 - 0.0001j;
ur = 1;
thk = 2.5;
f = 30;
k0 = f * (2*pi/299.792458);

%% nLayer
wgA = 7.112;
wgB = 3.556;

%% Get Spectrum Functions
[specX1, specY1] = nLayerOpenEnded.getSpectrumRectangular(...
    wgA, wgB, 1, 0, "TE");

%% Integral Old
tic;
[phi(1, 1, :), phi_weights(1, 1, :)] = fejer2(1000, 0, pi/2);
Ie = 4 * integral(@(x) gamE(x, k0, er, ur, thk) .* ...
    x .* sum(specY1(x.*cos(phi), x.*sin(phi)).^2 .* sin(phi).^2 .* phi_weights, 3), ...
    0, inf, RelTol=1e-6);
Ih = 4 * integral(@(x) gamH(x, k0, er, ur, thk) .* ...
    x .* sum(specY1(x.*cos(phi), x.*sin(phi)).^2 .* cos(phi).^2 .* phi_weights, 3), ...
    0, inf, RelTol=1e-6);
toc;

%% Integral 1
[krc1, krc_weights_E1, krc_weights_H1] = getContourWeights(...
    Nc, Nrho, Nphi, L, Lc, wgA, wgB);

I1e = sum(krc_weights_E1 .* gamE(krc1, k0, er, ur, thk));
I1h = sum(krc_weights_H1 .* gamH(krc1, k0, er, ur, thk));

%% Error
errE1 = db(I1e - Ie)
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
plot(1:numel(krc1), real(krc_weights_E1), "", LineWidth=1.5);
hold on;
plot(1:numel(krc1), imag(krc_weights_E1), "", LineWidth=1.5);grid on;
title("gamE Weights");

figure;
plot(1:numel(krc1), real(krc_weights_H1), "", LineWidth=1.5);
hold on;
plot(1:numel(krc1), imag(krc_weights_H1), "", LineWidth=1.5);grid on;
title("gamH Weights");





