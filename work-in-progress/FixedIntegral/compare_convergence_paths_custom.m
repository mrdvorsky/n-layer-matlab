clc;
clear;
close all;

%% Inputs
Nc = 11;
L = 5;
Lc = 4;

Lch = Lc;
Lcw = 40;


Nrho = 8000;
Nphi = 50;

%% Structure
er = {2 - 0.000001j};
ur = {1};
thk = {2.5};
f = 30;
k0 = f * (2*pi/299.792458);

%% nLayer
wgA = 7.112;
wgB = 3.556;

%% Contour Paths
% krToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + (x));
% krToKrcPrime = @(x) 1 + 1j .* Lc.^2 ./ ((x) + Lc).^2;

% krToKrc = @(x) x + 1j .* Lc .* x ./ sqrt(Lc + (x).^2);
% krToKrcPrime = @(x) 1 + 1j .* Lc.^2 ./ ((x).^2 + Lc).^1.5;
% 
% krToKrc = @(x) x + 1j*min(x, 1);
% krToKrcPrime = @(x) 1 + 1j*(x <= 1);

% krToKrc = @(x) x + 1j*tanh(x);
% krToKrcPrime = @(x) 1 + 1j*sech(x).^2;

% krToKrc = @(x) x + 1j .* x .* exp(-0.5 * (x ./ Lc).^2);
% krToKrcPrime = @(x) 1 + 1j .* exp(-0.5 * (x ./ Lc).^2) .* (1 - (x ./ Lc).^2);
% 
% krToKrc = @(x) x + 1j .* x .* exp(-0.5 * (x ./ Lc).^2);
% krToKrcPrime = @(x) 1 + 1j .* exp(-0.5 * (x ./ Lc).^2) .* (1 - (x ./ Lc).^2);


% krToKrcImag = @(x) Lch * tanh(Lch * x);
% krToKrcPrimeImag = @(x) Lch.^2 * sech(Lch * x).^2;

% krToKrcImag = @(x) Lch * x ./ sqrt(Lch + x.^2);
% krToKrcPrimeImag = @(x) Lch.^2 * x ./ (Lch + x.^2).^1.5;

% krToKrc = @(x) x + 1j * krToKrcImag(x);
% krToKrcPrime = @(x) 1 + 1j * krToKrcPrimeImag(x);

% krToKrc = @(x) x + 1j * (krToKrcImag(x) - 0.5*krToKrcImag(0.1*x - Lcw) - 0.5*Lch);
% krToKrcPrime = @(x) 1 + 1j * (krToKrcPrimeImag(x) - 0.05*krToKrcPrimeImag(0.1*x - Lcw));
% 
% krToKrc = @(x) x + 1j * Lch * tanh(Lch * x) ./ sqrt(1 + (x ./ Lcw).^2);
% krToKrcPrime = @(x) 1 + 1j * (Lch.^2 * (Lcw + x) .* sech(Lch * x).^2 - Lch * tanh(Lch * x)) ...
%     ./ (Lcw + x).^2;

krToKrc = @(x) x + 1j * (Lch.*Lcw) * x ./ sqrt((Lch + x.^2).*(3*Lcw.^2 + x.^2));
krToKrcPrime = @(x) 1 + 1j * (Lch.*Lcw) * (3*Lch*Lcw.^2 - x.^4) ./ ((Lch + x.^2).*(3*Lcw.^2 + x.^2)).^1.5;

% krToKrc = @(x) x + 1j .* Lc.^2 .* x ./ (Lc.^2 + x.^2);
% krToKrcPrime = @(x) 1 + 1j * Lc.^2 .* (Lc.^2 - x.^2) ./ (Lc.^2 + x.^2).^2;

% ezplot(@(x) imag(krToKrcPrime(x)), [0, 20], 20); ylim([-inf, inf]);
% ezplot(@(x) imag(krToKrc(x)), [0, 20], 21); ylim([-inf, inf]);
% return;

%% Integral Old
NL = nLayerRectangular(1, 0, waveguideBand="Ka");
specY1 = NL.modeStructs.EySpec;

[phi(1, 1, :), phi_weights(1, 1, :)] = fejer2(50, 0, pi/2);
we = @(x) gamE(x, k0, er, ur, thk) .* ...
    x .* sum(specY1(x.*cos(phi), x.*sin(phi)).^2 .* sin(phi).^2 .* phi_weights, 3);
wh = @(x) gamH(x, k0, er, ur, thk) .* ...
    x .* sum(specY1(x.*cos(phi), x.*sin(phi)).^2 .* cos(phi).^2 .* phi_weights, 3);

tic;
Ie = 4 * integral(we, 0, inf, RelTol=1e-8);
Ih = 4 * integral(wh, 0, inf, RelTol=1e-8);

I2e = 4 * integral(@(x) we(krToKrc(x)) .* krToKrcPrime(x), 0, inf, RelTol=1e-8);
I2h = 4 * integral(@(x) wh(krToKrc(x)) .* krToKrcPrime(x), 0, inf, RelTol=1e-8);
toc;

%% Integral 1
[krc1, krc_weights_E1, krc_weights_H1, krc_error_E1, krc_error_H1] = getContourWeights_customPath(...
    Nc, Nrho, Nphi, L, wgA, wgB, krToKrc, krToKrcPrime);

I1e = sum(krc_weights_E1 .* gamE(krc1, k0, er, ur, thk));
I1h = sum(krc_weights_H1 .* gamH(krc1, k0, er, ur, thk));

%% Integral 2
% [krc2, krc_weights_E2, krc_weights_H2] = getContourWeights(...
%     Nc, Nrho, Nphi, L, Lc, wgA, wgB);
% 
% I2e = sum(krc_weights_E2 .* gamE(krc2, k0, er, ur, thk));
% I2h = sum(krc_weights_H2 .* gamH(krc2, k0, er, ur, thk));

%% Error
errE1 = db(I1e - Ie)
errH1 = db(I1h - Ih)

% errE2 = db(I2e - Ie)
% errH2 = db(I2h - Ih)
% 
% errEest = db(sum(krc_error_E1 .* gamE(krc1, k0, er, ur, thk)))
% errHest = db(sum(krc_error_H1 .* gamE(krc1, k0, er, ur, thk)))

%% Integral Compare
% krInf = @(kr) L * cot(0.5 * (pi - kr*pi)).^2;
% gH = @(kr) gamH(krToKrc(kr), k0, er, ur, thk) ./ sqrt(1 + kr.^2);
% gE = @(kr) gamE(krToKrc(kr), k0, er, ur, thk) .* sqrt(1 + kr.^2);
% 
% gHfin = @(kr) gH(krInf(kr));
% gEfin = @(kr) gE(krInf(kr));
% 
% [x, w, w2] = fejer2(Nc, 0, 1);
% db(integral(gHfin, 0, 1) - sum(gHfin(x).*w))
% db(integral(gEfin, 0, 1) - sum(gEfin(x).*w))
% db(sum(gHfin(x).*w2))
% db(sum(gEfin(x).*w2))


%% Plot Gamma
[krPlot] = fejer2_halfOpen(10000, L);
krcPlot = krToKrc(krPlot);

figure;
plot(1:numel(krcPlot), real(gamH(krcPlot, k0, er, ur, thk) ./ sqrt(1 + krPlot.^2)), "", LineWidth=1.5);
hold on;
plot(1:numel(krcPlot), imag(gamH(krcPlot, k0, er, ur, thk) ./ sqrt(1 + krPlot.^2)), "", LineWidth=1.5);
grid on;
title("gamH");

figure;
plot(1:numel(krcPlot), real(gamE(krcPlot, k0, er, ur, thk) .* sqrt(1 + krPlot.^2)), "", LineWidth=1.5);
hold on;
plot(1:numel(krcPlot), imag(gamE(krcPlot, k0, er, ur, thk) .* sqrt(1 + krPlot.^2)), "", LineWidth=1.5);
grid on;
title("gamE");

figure;
plot(1:numel(krc1), real(krc_weights_E1), "", LineWidth=1.5);
hold on;
plot(1:numel(krc1), imag(krc_weights_E1), "", LineWidth=1.5);grid on;
title("gamE Weights");

figure;
plot(1:numel(krc1), real(krc_weights_H1), "", LineWidth=1.5);
hold on;
plot(1:numel(krc1), imag(krc_weights_H1), "", LineWidth=1.5);grid on;
title("gamH Weights");



figure;
plot(1:numel(krcPlot), real(we(krcPlot) ./ sqrt(1 + (krcPlot).^2)), "", LineWidth=1.5);
hold on;
plot(1:numel(krcPlot), real(we(krcPlot) ./ sqrt(1 + (krcPlot).^2)), "", LineWidth=1.5);grid on;
title("Mode Fun E");

figure;
plot(1:numel(krcPlot), real(wh(krcPlot) .* sqrt(1 + (krcPlot).^2)), "", LineWidth=1.5);
hold on;
plot(1:numel(krcPlot), real(wh(krcPlot) .* sqrt(1 + (krcPlot).^2)), "", LineWidth=1.5);grid on;
title("Mode Fun H");



