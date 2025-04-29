clc;
clear;
% close all;

%% Inputs
Lc = 0.1;

ccL = 0.3;
ccN = 5;

nX = 1000;
nY = 1000;

nR = 8000;

nS = 4000;

%% nLayer
NL = nLayerRectangular(1, 0, waveguideBand="x", convergenceAbsTol=1e-7);

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

%% Helper
% Nphi = 400;
% [phi1(1, 1, 1, :), phi1_weights(1, 1, 1, :), phi1_error(1, 1, 1, :)] = fejer2(Nphi, 0, 2*pi);

phi1(1, 1, 1, :) = (((1:Nphi) - 1) ./ Nphi) .* (2*pi);
phi1_weights(1, 1, 1, :) = 0*phi1 + (2*pi ./ Nphi);
phi1_error = phi1_weights;
phi1_error(1:2:end) = -phi1_error(1:2:end);

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
[phi2(1, :), phi2_weights(1, :)] = fejer2(400, 0, 2*pi);

r2(:, 1) = linspace(0, convR, nR);

xCirc = r2 .* cos(phi2);
yCirc = r2 .* sin(phi2);

tic;
EyyCirc = reshape(...
    nufftn(specEyy, {kx(:), ky(:)}, [xCirc(:), yCirc(:)] ./ (2*pi)), ...
    numel(r2), numel(phi2));
toc;

Eyy0 = sum(phi2_weights .* EyyCirc, 2);
Eyy2 = sum(cos(2*phi2) .* phi2_weights .* EyyCirc, 2);

%% Power Series 2
Eyy0_gi = griddedInterpolant(r2(1:2:end).^2, Eyy0(1:2:end), "spline");
Eyy0_fit = Eyy0_gi(r2.^2);

Eyy2_gi = griddedInterpolant(r2(1:2:end).^2, Eyy2(1:2:end), "spline");
Eyy2_fit = Eyy2_gi(r2.^2);

%% Plotting
figure;
showImage(xSamp, ySamp, EyySamp, DisplayFormat="Magnitude");

figure;
plot(r2, real(Eyy0), "", LineWidth=1.5);
hold on;
plot(r2, real(Eyy0_fit), ":", LineWidth=1.5);
grid on;

figure;
plot(r2, db((Eyy0 - Eyy0_fit) ./ max(abs(Eyy0))), ".", LineWidth=1.5);
grid on;
ylim([-300, 0]);


figure;
plot(r2, real(Eyy2), "", LineWidth=1.5);
hold on;
plot(r2, real(Eyy2_fit), ":", LineWidth=1.5);
grid on;

figure;
plot(r2, db((Eyy2 - Eyy2_fit) ./ max(abs(Eyy2))), ".", LineWidth=1.5);
grid on;
ylim([-300, 0]);




