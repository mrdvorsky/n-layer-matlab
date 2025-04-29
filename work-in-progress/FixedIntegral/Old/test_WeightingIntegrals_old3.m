clc;
clear;
close all;

%% Inputs
Lc = 0.1;

ccL = 0.3;
ccN = 5;

nX = 1000;
nY = 1000;

nR = 16000;

nS = 4000;

%% nLayer
NL = nLayerRectangular(1, 0, waveguideBand="x", convergenceAbsTol=1e-7);

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

%% Helper
[phi1(1, 1, 1, :), phi1_weights(1, 1, 1, :)] = fejer2(400, 0, 2*pi);

modeFunH_kr = @(kr) sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
    .* cos(phi1).^2 .* phi1_weights, 4);
modeFunH = @(kr) kr .* modeFunH_kr(kr);
modeFunE_kr = @(kr) sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
    .* sin(phi1).^2 .* phi1_weights, 4);
modeFunE = @(kr) kr .* modeFunE_kr(kr);

xToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + x);
xToKrc_weights = @(x) 1 + 1j .* Lc.^2 ./ (x + Lc).^2;

%% Compute Convolutions
convR = hypot(NL.waveguideA, NL.waveguideB);
xSamp = linspace(-convR, convR, nX);
ySamp = linspace(-convR, convR, nY);

[kx, ky] = fftCoordinates(xSamp, ySamp, ApplyFftShift=true);
specEyy = specFunEy(kx, ky).^2;

scaleFactor = numel(specEyy) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ./ (2*pi);
EyySamp = scaleFactor .* fftshift(ifft2(ifftshift(specEyy)));

%% Nonuniform Sampling
[phi2(1, :), phi2_weights(1, :)] = fejer2(100, 0, 2*pi);
[phir2(:, 1), phir2_weights(:, 1)] = fejer2(nR, 0, pi);
% phir2 = [phir2; phir2 + pi];
% phir2_weights = [phir2_weights; -phir2_weights];

r2(:, 1) = sin(phir2) .* convR;

xCirc = r2 .* cos(phi2);
yCirc = r2 .* sin(phi2);

tic;
EyyCirc = reshape(...
    nufftn(specEyy, {kx(:), ky(:)}, [xCirc(:), yCirc(:)] ./ (2*pi)), ...
    numel(r2), numel(phi2));
toc;

Eyy0 = sum(phi2_weights .* EyyCirc, 2);
Eyy2 = sum(cos(2*phi2) .* phi2_weights .* EyyCirc, 2);

%% Power Series
% W0 = (1 - (r2./convR).^2).^(0:nS - 1);
% cS0 = W0 \ Eyy0;
% 
% W2 = (r2 ./ convR).^2 .* (1 - (r2./convR).^2).^(0:nS - 1);
% cS2 = W2 \ Eyy2;
% 
% Eyy0_fit = W0 * cS0;
% Eyy2_fit = W2 * cS2;

%% Power Series 2
k(1, :) = 2*(0:nS - 1);

gm0 = 2 * sum(cos(k .* phir2) .* phir2_weights .* Eyy0, 1) ./ pi;
gm0(k == 0) = 0.5 * gm0(k == 0);

Eyy0_fit = sum(cos(k .* phir2) .* gm0, 2);

gm2 = 2 * sum(cos(k .* phir2) .* phir2_weights .* Eyy2 ./ sin(phir2).^2, 1) ./ pi;
gm2(k == 0) = 0.5 * gm2(k == 0);

Eyy2_fit = sum(cos(k .* phir2) .* gm2, 2) .* sin(phir2).^2;

%% Plotting
figure;
showImage(xSamp, ySamp, EyySamp, DisplayFormat="Magnitude");

figure;
plot(phir2, real(Eyy0), "", LineWidth=1.5);
hold on;
plot(phir2, real(Eyy0_fit), ":", LineWidth=1.5);
grid on;

figure;
plot(phir2, db(Eyy0 - Eyy0_fit), "", LineWidth=1.5);
grid on;


figure;
plot(phir2, real(Eyy2) ./ sin(phir2).^2, "", LineWidth=1.5);
hold on;
plot(phir2, real(Eyy2_fit), ":", LineWidth=1.5);
grid on;

figure;
plot(phir2, db(Eyy2 - Eyy2_fit), "", LineWidth=1.5);
grid on;





