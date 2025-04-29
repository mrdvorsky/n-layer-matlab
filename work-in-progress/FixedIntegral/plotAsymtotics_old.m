clc;
clear;
close all;

%% Inputs
krMax = 1000;
Nkr = 20000;

Nphi = 800;

%% nLayer
NL.waveguideA = 1;
NL.waveguideB = 1/2;

%% Get Spectrum Functions
[spEx1, spEy1] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");
[spEx2, spEy2] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 3, 0, "TE");

[phi(1, 1, :), phi_weights(1, 1, :)] = fejer2(Nphi, 0, 0.5*pi);
phi_weights = 4*phi_weights;

modeFunE = @(kr) kr .* sum(phi_weights ...
    .* (sin(phi) .* spEx1(kr .* cos(phi), kr .* sin(phi)) ...
    - cos(phi) .* spEy1(kr .* cos(phi), kr .* sin(phi))) ...
    .* (sin(phi) .* spEx2(kr .* cos(phi), kr .* sin(phi)) ...
    - cos(phi) .* spEy2(kr .* cos(phi), kr .* sin(phi))), ...
    3);
modeFunH = @(kr) kr .* sum(phi_weights ...
    .* (cos(phi) .* spEx1(kr .* cos(phi), kr .* sin(phi)) ...
    + sin(phi) .* spEy1(kr .* cos(phi), kr .* sin(phi))) ...
    .* (cos(phi) .* spEx2(kr .* cos(phi), kr .* sin(phi)) ...
    + sin(phi) .* spEy2(kr .* cos(phi), kr .* sin(phi))), ...
    3);
modeFunEH = @(kr) modeFunE(kr) + modeFunH(kr);

%% Sample
krq(:, 1) = linspace(0, krMax, Nkr);
specq1 = modeFunEH(krq);

specq2 = besselj(0.5, 0.5 ./ (2) .* krq).^2 ./ krq;

%% Plotting
figure;
semilogx(krq, db(specq1), "", LineWidth=1.5);
grid on;
title("Spectrum 1");
yLims = ylim();

figure;
semilogx(krq, db(specq2), "", LineWidth=1.5);
grid on;
title("Spectrum 2");
ylim(yLims);

figure;
semilogx(krq, db(specq1 - specq2), "", LineWidth=1.5);
grid on;
title("Spectrum Error");


figure;
plot(krq, real(specq1), "", LineWidth=1.5);
hold on;
plot(krq, specq2, "", LineWidth=1.5);
grid on;
title("Spectrum 1");
xlim(krMax * [0.99, 1]);




