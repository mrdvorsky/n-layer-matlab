clc;
clear;
close all;

%% Inputs
Lc = 1.5;
ccL = 1.1;

ccN = 20;

N1 = 50;
N2 = 40;

er = 2 - 0.1j;
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
alpha = @(x) x + 1j .* Lc .* x ./ (Lc + abs(x));
alphaP = @(x) 1 + 1j .* Lc.^2 ./ (abs(x) + Lc).^2;

beta = @(x) 0.5 * (sqrt(x.^2 + (2-2j)*x + 2j) + x - 1 - 1j);
betaP = @(x) 0.5 * (1 + (x + 1 - 1j) ...
    ./ sqrt(x.^2 + (2-2j)*x + 2j));

%% Helper 2
phi1(1, 1, 1, :) = 2*pi * ((1:4*N1) - 0.5) ./ (4*N1);
phi1_weights = 0*phi1(1:N1) + 4*(2*pi ./ (4*N1));
phi1 = phi1(1:N1);

% phi2(1, :) = 2*pi * ((1:4*N2) - 0.5) ./ (4*N2);
% phi2_weights = 0*phi2 + 1*(2*pi ./ (4*N2));
% % phi2 = phi2(1:N2);

% [phi2(1, :), phi2_weights(1, :)] = fejer2(N2, 0, pi/2);
% phi2_weights = 4*phi2_weights;

modeFunE = @(kr) kr .* sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
    .* sin(phi1).^2 .* phi1_weights, 4);
modeFunH = @(kr) kr .* sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
    .* cos(phi1).^2 .* phi1_weights, 4);

%% Compute Weighting 1
% tic;
% [x(1, :), kr_weights_E(1, :)] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
%     WeightingFunction=@(x) modeFunE(alpha(x)) .* alphaP(x));
% kr = alpha(x);
% toc;
% 
% tic;
% [~, kr_weights_H(1, :)] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
%     WeightingFunction=@(x) modeFunH(alpha(x)) .* alphaP(x));
% toc;

tic;
[x(1, :), kr_weights_E(1, :)] = fejer2_halfOpen(ccN, ccL, ...
    WeightingFunction=@(x) modeFunE(alpha(x)) .* alphaP(x) .* (1 + 0*sqrt(1 + alpha(x).^2)));
kr = alpha(x);
% kr_weights_E = kr_weights_E ./ sqrt(1 + alpha(x).^2);
toc;

tic;
[~, kr_weights_H(1, :)] = fejer2_halfOpen(ccN, ccL, ...
    WeightingFunction=@(x) modeFunH(alpha(x)) .* alphaP(x) .* sqrt(1 + alpha(x).^2));
kr_weights_H = kr_weights_H ./ sqrt(1 + alpha(x).^2);
toc;

%% Compare Integral
I1e = integral(@(x) modeFunE(x) .* gamE(x, k0, er, ur, thk), 0, inf);
I1h = integral(@(x) modeFunH(x) .* gamH(x, k0, er, ur, thk), 0, inf);

I2e = sum(kr_weights_E .* gamE(kr, k0, er, ur, thk));
I2h = sum(kr_weights_H .* gamH(kr, k0, er, ur, thk));


errE2 = db(I1e - I2e)
errH2 = db(I1h - I2h)


%% Plot
figure;
plot(1:numel(kr), real(gamH(kr, k0, er, ur, thk) ./ sqrt(1 + kr.^2)), "", LineWidth=1.5);
hold on;
plot(1:numel(kr), imag(gamH(kr, k0, er, ur, thk) ./ sqrt(1 + kr.^2)), "", LineWidth=1.5);
grid on;
legend("Real", "Imag");
title("H");

figure;
plot(1:numel(kr), real(gamE(kr, k0, er, ur, thk) ./ (1 + 0*sqrt(1 + kr.^2))), "", LineWidth=1.5);
hold on;
plot(1:numel(kr), imag(gamE(kr, k0, er, ur, thk) ./ (1 + 0*sqrt(1 + kr.^2))), "", LineWidth=1.5);
grid on;
legend("Real", "Imag");
title("E");



