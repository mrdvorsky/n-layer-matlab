clc;
clear;
close all;

%% Inputs
Nc = 1;
L = 1;
Lc = 1.0;

Nrho = 16000;
Nphi = 950;

krMax = 500;

%% nLayer
NL.waveguideA = 7.112;
NL.waveguideB = 7.112/2;

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

[phi(1, 1, :), phi_weights(1, 1, :)] = fejer2(Nphi, 0, pi/2);
phi_weights = 4*phi_weights;

modeFunE = @(kr) (kr) .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* sin(phi).^2 .* phi_weights, 3);
modeFunH = @(kr) (kr) .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* cos(phi).^2 .* phi_weights, 3);
modeFunEH = @(kr) modeFunE(kr) + modeFunH(kr);

%% Sample
krQuery = linspace(-krMax, krMax, 2*Nrho - 1);

specEH = modeFunEH(abs(krQuery)) .* (1 + 1 * krQuery.^2).^1;
rQuery = fftCoordinates(krQuery, ApplyFftShift=true);

scaleFactor = abs(krQuery(2) - krQuery(1)) * numel(krQuery);
spatEH = fftshift(ifft(specEH)) * scaleFactor;
spatEH = spatEH ./ max(abs(spatEH));

%% Plotting
figure;
plot(krQuery, db(specEH), "", LineWidth=1.5);
grid on;
title("Spectrum");

figure;
plot(rQuery, db(spatEH), "", LineWidth=1.5);
grid on;
title("Spatial");






