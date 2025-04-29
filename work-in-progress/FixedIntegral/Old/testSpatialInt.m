clc;
clear;
close all;

%% Inputs
krMax = 400;
Nkr = 32000;

Nphi = 800;

%% nLayer
NL.waveguideA = 7.112;
NL.waveguideB = 7.112/2;

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

[phi(1, 1, :), phi_weights(1, 1, :)] = fejer2(Nphi, 0, pi/2);
phi_weights = 4*phi_weights;

modeFunE = @(kr) kr .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* sin(phi).^2 .* phi_weights, 3);
modeFunH = @(kr) kr .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* cos(phi).^2 .* phi_weights, 3);
modeFunEH = @(kr) modeFunE(kr) + modeFunH(kr);

%% Test Integral
[x, w] = fejer2_halfOpen(Nkr, 1);
I1 = sum(modeFunEH(x) .* w);

%% Sample
krQuery(:, 1) = linspace(-krMax, krMax, 2*Nkr - 1);
specQuery = modeFunEH(abs(krQuery));

rQuery = fftCoordinates(krQuery, ApplyFftShift=true);
scaleFactor = abs(krQuery(2) - krQuery(1)) .* numel(krQuery);
spatQuery = fftshift(ifft(ifftshift(specQuery))) .* scaleFactor;

%% Plotting
figure;
semilogx(krQuery, db(specQuery), "", LineWidth=1.5);
grid on;
title("Spectrum 1");

figure;
plot(rQuery, abs(spatQuery), "", LineWidth=1.5);
grid on;
title("Spatial");






