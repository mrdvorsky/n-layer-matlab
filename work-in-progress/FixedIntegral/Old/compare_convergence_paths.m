clc;
clear;
close all;

%% Inputs
Nc = 1;
L = 1;
Lc = 1.0;

Nrho = 16000;
Nphi = 2*950;

%% Structure
er = 2 - 0.0001j;
ur = 1;
thk = 2.5;
f = 30;
k0 = f * (2*pi/299.792458);

%% nLayer
% NL = nLayerRectangular(1, 0, waveguideBand="ka", convergenceAbsTol=1e-7);
NL.waveguideA = 7.112;
NL.waveguideB = 7.112/2;

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");


%% Integral Accurate
[krc, krc_weights_E1, krc_weights_H1] = getContourWeights(...
    Nc, Nrho, Nphi, L, Lc, NL.waveguideA, NL.waveguideB);

I1e = sum(krc_weights_E1 .* gamE(krc, k0, er, ur, thk));
I1h = sum(krc_weights_H1 .* gamH(krc, k0, er, ur, thk));

%% Integral 2
[krc, krc_weights_E2, krc_weights_H2] = getContourWeights(...
    Nc, 1*Nrho, 4*Nphi, L, Lc, NL.waveguideA, NL.waveguideB);

I2e = sum(krc_weights_E2 .* gamE(krc, k0, er, ur, thk));
I2h = sum(krc_weights_H2 .* gamH(krc, k0, er, ur, thk));

%% Error
errEH2 = db(I1e + I2h - I2e - I2h)


errE2 = db(I1e - I2e)
errH2 = db(I1h - I2h)


errSum_db1 = db(1 - sum(krc_weights_E1 + krc_weights_H1))
errSum_db2 = db(1 - sum(krc_weights_E2 + krc_weights_H2))








