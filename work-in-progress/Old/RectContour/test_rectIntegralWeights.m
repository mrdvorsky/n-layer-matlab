clc;
clear;
close all;

%% Inputs
f(:, 1) = linspace(8.2, 12.4, 21);

er = [2 - 0.001j];
ur = [1];
thk = [10];

Lc = 0.1;

ccOrder = 50;
ccL = 0.3;

%% nLayer
NL = nLayerRectangular(1, 0, waveguideBand="x", convergenceAbsTol=1e-7);

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

%% Calculate
gam1 = NL.calculate(f, er, ur, thk);

%% Helper
k0 = 2*pi* f ./ NL.speedOfLight;

[phi(1, 1, 1, :), phi_weights(1, 1, 1, :)] = fejer2(50, 0, pi/2);

modeFunH = @(kr) 4*kr .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* cos(phi).^2 .* phi_weights, 4);
modeFunE = @(kr) 4*kr .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* sin(phi).^2 .* phi_weights, 4);

xToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + x);
xToKrc_weights = @(x) 1 + 1j .* Lc.^2 ./ (x + Lc).^2;

%% Integral Complex Weighted
tic;
[x, x_weights_H] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
    WeightingFunction=@(x) modeFunH(xToKrc(x)) .* xToKrc_weights(x) .* (1 + x));
toc;
[~, x_weights_E] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
    WeightingFunction=@(x) modeFunE(xToKrc(x)) .* xToKrc_weights(x));
toc;

for ff = 1:numel(k0)
    gamE = @(kr) out2(@(k) nLayerOpenEnded.computeGamma0(k, k0(ff), er, ur, thk), kr);
    gamH = @(kr) nLayerOpenEnded.computeGamma0(kr, k0(ff), er, ur, thk);

    valH(ff, 1) = sum(gamH(xToKrc(x)) .* x_weights_H ./ (1 + x));
    valE(ff, 1) = sum(gamE(xToKrc(x)) .* x_weights_E);
end

A = valH + valE;
K = k0 ./ sqrt(k0.^2 - NL.cutoffBeta_TE.^2);

gam4 = (1 - A.*K) ./ (1 + A.*K);

%% Error
err4_db = db(max(gam1 - gam4))

%% Plotting
figure;
plot(gam1, "", LineWidth=1.5);
hold on;
% plot(gam2, "o", LineWidth=1.5);
% plot(gam3, "x", LineWidth=1.5);
plot(gam4, ".", LineWidth=1.5);
zplane([]);




