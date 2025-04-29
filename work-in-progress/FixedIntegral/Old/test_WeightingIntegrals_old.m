clc;
clear;
close all;

%% Inputs
Lc = 0.1;

ccL = 0.3;
ccN = 5;

xs(:, 1) = linspace(-160, 160, 2001);
ys(1, :) = linspace(-160, 160, 2001);

rs(:, 1) = linspace(0.001, 40, 4001);

kr(:, 1) = linspace(0.001, 5, 1001);

%% nLayer
NL = nLayerRectangular(1, 0, waveguideBand="x", convergenceAbsTol=1e-7);

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

%% Helper
[phi1(1, 1, 1, :), phi1_weights(1, 1, 1, :)] = fejer2(400, 0, 2*pi);

modeFunH = @(kr) kr .* sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
    .* cos(phi1).^2 .* phi1_weights, 4);
% modeFunE = @(kr) kr .* sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
%     .* sin(phi1).^2 .* phi1_weights, 4);

% xToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + x);
% xToKrc_weights = @(x) 1 + 1j .* Lc.^2 ./ (x + Lc).^2;

%% Compute Weighting
% [x, x_weights_E] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
%     WeightingFunction=@(x) modeFunE(xToKrc(x)) .* xToKrc_weights(x));

% for ii = 1:ccN
%     wf = @(x) modeFunE(xToKrc(x)) .* xToKrc_weights(x);
% 
%     moments1(ii) = 0.5 ./ ccL .* integral(...
%         @(x) wf(x) ...
%         .* cos(2*(ii - 1) * acot(sqrt(x ./ ccL))), ...
%         0, inf);
% end

%% Spatial
[kx, ky] = fftCoordinates(xs, ys);

specEy = specFunEy(kx, ky) + 0*kx.*ky;
specConv = specEy.^2 .* cos(angle(kx + 1j.*ky)).^2;

scaleFactor = numel(specEy) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
        ./ (2*pi);
gamH = scaleFactor .* ifftshift(ifft2(specConv .* (hypot(kx, ky) < max(kx))));
gamH_interp = griddedInterpolant({xs, ys}, gamH, "spline");

[phi2(1, :), phi2_weights(1, :)] = fejer2(400, 0, 2*pi);
gamH_r2 = sum(gamH_interp(rs .* cos(phi2), rs .* sin(phi2)) .* phi2_weights, 2);

%% Compute GamH_r2
krSamp(:, 1) = linspace(0, 40, 10001);
gamH_r1 = sum(modeFunH(krSamp).' .* besselj(0, rs .* krSamp.'), 2) ...
    .* abs(krSamp(2) - krSamp(1));

%% Compute Spectrums to Compare
spec1 = modeFunH(kr);
spec2 = sum(gamH_r1.' .* besselj(0, rs.' .* kr) .* rs.', 2) ...
     .* abs(rs(2) - rs(1)) .* kr;

%% Plotting
figure;
plot(kr, abs(spec1), "", LineWidth=1.5);
grid on;

figure;
plot(kr, abs(spec2), "", LineWidth=1.5);
grid on;

% figure;
% plot(rs, abs(gamH_r1), "", LineWidth=1.5);
% grid on;
% 
% figure;
% plot(rs, abs(gamH_r2), "", LineWidth=1.5);
% grid on;

% figure;
% showImage(xs, ys, gamH, DisplayFormat="Magnitude");

%% Old
% tic;
% [x, x_weights_H] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
%     WeightingFunction=@(x) modeFunH(xToKrc(x)) .* xToKrc_weights(x) .* (1 + x));
% toc;
% [~, x_weights_E] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
%     WeightingFunction=@(x) modeFunE(xToKrc(x)) .* xToKrc_weights(x));
% toc;
% 
% for ff = 1:numel(k0)
%     gamE = @(kr) out2(@(k) nLayerOpenEnded.computeGamma0(k, k0(ff), er, ur, thk), kr);
%     gamH = @(kr) nLayerOpenEnded.computeGamma0(kr, k0(ff), er, ur, thk);
% 
%     valH(ff, 1) = sum(gamH(xToKrc(x)) .* x_weights_H ./ (1 + x));
%     valE(ff, 1) = sum(gamE(xToKrc(x)) .* x_weights_E);
% end






