clc;
clear;
close all;

%% Inputs
Nc = 50;
L = 1;
Lc = 1.0;

Nrho = 16000;
Nphi = 950;

%% Structure
er = 2 - 0.0001j;
ur = 1;
thk = 2.5;
f = 30;
k0 = f * (2*pi/299.792458);

%% nLayer
NL = nLayerRectangular(1, 0, waveguideBand="ka", convergenceAbsTol=1e-7);

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

%% Helper
[phi(1, 1, :), phi_weights(1, 1, :)] = fejer2(Nphi, 0, pi/2);
phi_weights = 4*phi_weights;

modeFunE = @(kr) kr .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* sin(phi).^2 .* phi_weights, 3);
modeFunH = @(kr) kr .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* cos(phi).^2 .* phi_weights, 3);

%% Compare Integral
I1e = integral(@(x) modeFunE(x) .* gamE(x, k0, er, ur, thk), 0, inf);
I1h = integral(@(x) modeFunH(x) .* gamH(x, k0, er, ur, thk), 0, inf);

[krc, krc_weights_E, krc_weights_H] = getContourWeights(...
    Nc, Nrho, Nphi, L, Lc, NL.waveguideA, NL.waveguideB);
I2e = sum(krc_weights_E .* gamE(krc, k0, er, ur, thk));
I2h = sum(krc_weights_H .* gamH(krc, k0, er, ur, thk));

errEH2 = db(I1e + I2h - I2e - I2h)


errE2 = db(I1e - I2e)
errH2 = db(I1h - I2h)

errSum_db2 = db(1 - sum(krc_weights_E + krc_weights_H))








