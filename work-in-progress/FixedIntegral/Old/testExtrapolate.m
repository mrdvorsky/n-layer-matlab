clc;
clear;
close all;

%% Inputs
krMax1 = 2000;
Nkr1 = 16000;

krMax2 = 200;
Nkr2 = 6400;
Lt = 1;

Nphi = 800;

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
krQuery1(:, 1) = linspace(0, krMax1, Nkr1);
specQuery1 = modeFunEH(krQuery1);

specQuery2 = extrapolateSpectrum(modeFunEH, Lt, krMax2, Nkr2, krQuery1);

[x, w] = fejer2_halfOpen(32000, 1);
Isum1 = db(1 - sum(extrapolateSpectrum(modeFunEH, Lt, krMax2, Nkr2, x) .* w))
Isum2 = db(1 - sum(modeFunEH(x) .* w))
Isum3 = db(1 - sum(specQuery1) * abs(krQuery1(2) - krQuery1(1)))

%% Plotting
figure;
semilogx(krQuery1, db(specQuery1), "", LineWidth=1.5);
grid on;
title("Spectrum 1");
yLims = ylim();

figure;
semilogx(krQuery1, db(specQuery2), "", LineWidth=1.5);
grid on;
title("Spectrum 2");
ylim(yLims);


figure;
semilogx(krQuery1, db(specQuery2 - specQuery1), "", LineWidth=1.5);
grid on;
title("Error");

% figure;
% plot(rQuery, abs(spatEH), "", LineWidth=1.5);
% grid on;
% title("Spatial");






