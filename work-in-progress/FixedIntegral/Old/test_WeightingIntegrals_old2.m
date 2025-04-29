clc;
clear;
close all;

%% Inputs
Lc = 0.1;

ccL = 0.3;
ccN = 5;

xs(:, 1) = 320 * linspace(-1, 1, 4001);
ys(1, :) = 320 * linspace(-1, 1, 4001);

rs(:, 1) = 40 * linspace(0.0000001, 1, 2001);

kr(:, 1) = linspace(0.001, 5, 4001);

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
% modeFunE = @(kr) kr .* sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
%     .* sin(phi1).^2 .* phi1_weights, 4);

xToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + x);
xToKrc_weights = @(x) 1 + 1j .* Lc.^2 ./ (x + Lc).^2;

%% Spatial
[kx, ky] = fftCoordinates(xs, ys);

specEy = specFunEy(kx, ky) + 0*kx.*ky;
specConv = specEy.^2 .* cos(angle(kx + 1j.*ky)).^2;

scaleFactor = numel(specEy) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
        ./ (2*pi);
gamH = scaleFactor .* ifftshift(ifft2(specConv .* (hypot(kx, ky) < 1*max(kx))));
gamH_interp = griddedInterpolant({xs, ys}, gamH, "linear");

[phi2(1, :), phi2_weights(1, :)] = fejer2(400, 0, 2*pi);
gamH_r2 = sum(gamH_interp(rs .* cos(phi2), rs .* sin(phi2)) .* phi2_weights, 2);

%% Compute Spectrums to Compare
spec1 = modeFunH(kr);


sumJ0_r(1, 1, :) = rs;
sumJ0_weights(1, 1, :) = gamH_r2 .* rs;

modeFunH2 = @(kr) sum(sumJ0_weights .* besselj(0, sumJ0_r .* kr), 3) ...
     .* abs(sumJ0_r(2) - sumJ0_r(1)) .* kr;
spec2 = modeFunH2(kr);

modeFunH3_1 = @(kr) 0.5 * sum(sumJ0_weights .* besselh(0, 1, sumJ0_r .* kr), 3) ...
     .* abs(sumJ0_r(2) - sumJ0_r(1)) .* kr;
modeFunH3_2 = @(kr) 0.5 * sum(sumJ0_weights .* besselh(0, 2, sumJ0_r .* kr), 3) ...
     .* abs(sumJ0_r(2) - sumJ0_r(1)) .* kr;
spec3 = modeFunH3_1(kr) + modeFunH3_2(kr);

%% Compute Weighting
% [x, x_weights_E] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
%     WeightingFunction=@(x) modeFunE(xToKrc(x)) .* xToKrc_weights(x));

for ii = 1:ccN
    ii

    wf = @(x) modeFunH(xToKrc(x)) .* xToKrc_weights(x);
    moments1(ii) = 0.5 ./ ccL .* integral(...
        @(x) wf(x) ...
        .* cos(2*(ii - 1) * acot(sqrt(x ./ ccL))), ...
        0, inf);

    wf = @(x) modeFunH2(xToKrc(x)) .* xToKrc_weights(x);
    moments2(ii) = 0.5 ./ ccL .* integral(...
        @(x) wf(x) ...
        .* cos(2*(ii - 1) * acot(sqrt(x ./ ccL))), ...
        0, 10);

    wf = @(x) modeFunH3_1(xToKrc(x)) .* xToKrc_weights(x);
    moments3_1(ii) = 0.5 ./ ccL .* integral(...
        @(x) wf(x) ...
        .* cos(2*(ii - 1) * acot(sqrt(x ./ ccL))), ...
        0, 10);

    wf = @(x) modeFunH3_2(xToKrc(x)) .* xToKrc_weights(x);
    moments3_2(ii) = 0.5 ./ ccL .* integral(...
        @(x) wf(x) ...
        .* cos(2*(ii - 1) * acot(sqrt(x ./ ccL))), ...
        0, -10j);

    moments3 = moments3_1 + moments3_2;
end

%% Plotting
figure;
plot(kr, real(spec1), "", LineWidth=1.5);
hold on;
plot(kr, real(spec2), ":", LineWidth=1.5);
plot(kr, real(spec3), ".", LineWidth=1.5);
grid on;

figure;
plot(kr, db(spec1 - spec2), "", LineWidth=1.5);
hold on;
plot(kr, db(spec1 - spec3), "", LineWidth=1.5);
grid on;







