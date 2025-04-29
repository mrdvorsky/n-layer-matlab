clc;
clear;
close all;

%% Inputs
Lc = 0.5;
ccL = 2.0;

ccN = 50;

N1 = 50;
N2 = 40;

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
% alpha = @(x) x + 1j .* Lc .* x ./ (Lc + (x));
% alphaP = @(x) 1 + 1j .* Lc.^2 ./ ((x) + Lc).^2;
% 
% beta = @(x) 0.5 * (sqrt(x.^2 + (2-2j)*x + 2j) + x - 1 - 1j);
% betaP = @(x) 0.5 * (1 + (x + 1 - 1j) ...
%     ./ sqrt(x.^2 + (2-2j)*x + 2j));

beta = @(x) x - 1j .* Lc .* x ./ (Lc + (x));
betaP = @(x) 1 - 1j .* Lc.^2 ./ ((x) + Lc).^2;

alpha = @(x) 0.5 * (sqrt(x.^2 + Lc*(2+2j)*x - 2j*Lc.^2) + x - Lc + 1j*Lc);
alphaP = @(x) 0.5 * (1 + (x + Lc + 1j*Lc) ...
    ./ sqrt(x.^2 + Lc*(2+2j)*x - 2j*Lc.^2));

%% Helper 2
phi1(1, 1, 1, :) = 2*pi * ((1:4*N1) - 0.5) ./ (4*N1);
phi1_weights = 0*phi1(1:N1) + 4*(2*pi ./ (4*N1));
phi1 = phi1(1:N1);

modeFunE = @(kr) kr .* sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
    .* sin(phi1).^2 .* phi1_weights, 4);
modeFunH = @(kr) kr .* sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
    .* cos(phi1).^2 .* phi1_weights, 4);

%% Compute Weighting 1
wfunE = @(x) modeFunE(alpha(x)) .* alphaP(x);
wfunH = @(x) modeFunH(alpha(x)) .* alphaP(x) .* sqrt(1 + (x).^2);
tic;
for kk = 1:ccN
    % MkE(kk) = integral(@(x) (2*ccL) * wfunE(ccL * cot(0.5*x).^2) ...
    %                         .* sin(kk*x) ./ (1 - cos(x)).^2, 0, pi);
    % MkH(kk) = integral(@(x) (2*ccL) * wfunH(ccL * cot(0.5*x).^2) ...
    %                         .* sin(kk*x) ./ (1 - cos(x)).^2, 0, pi);

    wfunE2 = @(x) (ccL) * wfunE(ccL * x) ...
                            .* sin(2*kk .* acot(sqrt(x))) ...
                            ./ sin(2 .* acot(sqrt(x)));
    wfunH2 = @(x) (ccL) * wfunH(ccL * x) ...
                            .* sin(2*kk .* acot(sqrt(x))) ...
                            ./ sin(2 .* acot(sqrt(x)));

    % MkE(kk) = integral(@(x) wfunE2(x), 0, inf);
    % MkH(kk) = integral(@(x) wfunH2(x), 0, inf);

    % MkE(kk) = integral(@(x) wfunE2(beta(x)) .* betaP(x), 0, inf);
    % MkH(kk) = integral(@(x) wfunH2(beta(x)) .* betaP(x), 0, inf);

    MkE(kk) = integral(@(x) modeFunE(x) .* alphaP(beta(x)) ...
                            .* sin(2*kk .* acot(sqrt(beta(x) ./ ccL))) ...
                            ./ sin(2 .* acot(sqrt(beta(x) ./ ccL))) ...
                            .* betaP(x), 0, inf);
    MkH(kk) = integral(@(x) (ccL) * wfunH(ccL * beta(x)) ...
                            .* sin(2*kk .* acot(sqrt(beta(x)))) ...
                            ./ sin(2 .* acot(sqrt(beta(x)))) ...
                            .* betaP(x), 0, inf);
end
toc;

tic;
[x(1, :), kr_weights_E(1, :)] = fejer2_halfOpen(ccN, ccL, ...
    WeightingMoments=MkE);
kr = alpha(x);
toc;

tic;
[~, kr_weights_H(1, :)] = fejer2_halfOpen(ccN, ccL, ...
    WeightingMoments=MkH);
kr_weights_H = kr_weights_H ./ sqrt(1 + (x).^2);
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



