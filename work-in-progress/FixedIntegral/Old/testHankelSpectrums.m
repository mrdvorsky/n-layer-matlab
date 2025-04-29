clc;
clear;
close all;

%% Inputs
N = 5000;
krMax = 500;
Nphi = 1000;

Nj = 1000;
Lj = 200;
vj = 1;

%% nLayer
NL.waveguideA = 7.112;
NL.waveguideB = 7.112/2;

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

[phi1(1, 1, :), w(1, 1, :)] = fejer2(Nphi, 0, pi/2);
w = 4*w;
modeFun = @(kr) kr .* sum(cos(2*phi1) .* w .* specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2, 3);

%% Sample
jvm(1, :) = besselj_zerosVec(vj, Nj);
krH(1, :) = Lj .* jvm ./ jvm(end);

Ht = besselj(vj, jvm .* jvm.' ./ jvm(end)) ./ besselj(vj + 1, jvm).^2;
um = (2 * Lj.^2 ./ jvm(end).^2) .* sum(Ht .* modeFun(krH), 2);

kr = linspace(0, krMax, N);
spec = modeFun(kr);

%% Plotting
figure;
semilogx(kr, db(spec), "", LineWidth=1.5);
grid on;
title("Spectrum 1");

figure;
semilogx(kr, db(spec), "", LineWidth=1.5);
grid on;
title("Spectrum 1");








