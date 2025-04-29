clc;
clear;
close all;

%% Inputs
Lc = 0.1;
ccL = 0.3;

ccOrder = 10;

N1 = 400;
N2 = 100;

krComp1(:, 1) = linspace(-20, 20, 1601);
krComp2(:, 1) = linspace(-100, 100, 8001);

%% nLayer
NL = nLayerRectangular(1, 0, waveguideBand="x", convergenceAbsTol=1e-7);

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

%% Helper
phi1(1, 1, 1, :) = 2*pi * ((1:4*N1) - 0.5) ./ (4*N1);
phi1_weights = 0*phi1(1:N1) + 4*(2*pi ./ (4*N1));
phi1 = phi1(1:N1);

% modeFunH_kr = @(kr, ph) sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
%     .* cos(phi1).^2 .* phi1_weights, 4);
% modeFunH = @(kr) kr .* modeFunH_kr(kr);
modeFunE_kr = @(kr) sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
    .* sin(phi1).^2 .* phi1_weights, 4);
% modeFunE = @(kr) kr .* modeFunE_kr(kr);

% xToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + x);
% xToKrc_weights = @(x) 1 + 1j .* Lc.^2 ./ (x + Lc).^2;

xToKrc = @(x) x + 1j .* Lc .* x ./ sqrt(Lc.^2 + x.^2);
xToKrc_weights = @(x) 1 + 1j .* Lc.^3 ./ (x.^2 + Lc.^2).^(1.5);

% xToKrc = @(x) 0.5 * (sqrt(x.^2 + (2+2j)*Lc*x - 2j*Lc.^2) + x - Lc + 1j*Lc);
% xToKrc_weights = @(x) 0.5 * (1 + (x + Lc + 1j*Lc) ...
%     ./ sqrt(x.^2 + (2+2j)*Lc*x - 2j*Lc.^2));

%% Helper 2
% phi2(1, :) = 2*pi * ((1:4*N2) - 0.5) ./ (4*N2);
% phi2_weights = 0*phi2 + 1*(2*pi ./ (4*N2));
% % phi2 = phi2(1:N2);

[phi2(1, :), phi2_weights(1, :)] = fejer2(N2, 0, pi/2);
phi2_weights = 4*phi2_weights;

modeFunE_krphi = @(kr, ph) specFunEy(kr .* cos(ph), kr .* sin(ph)).^2 ...
    .* sin(ph).^2;

%% Compute Weighting 1
[x, x_weights_E1] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
    WeightingFunction=@(x) modeFunE_kr(xToKrc(x)) .* xToKrc_weights(x));

%% Compute Weighting 2
for pp = 1:numel(phi2)
    w = @(x) modeFunE_krphi(xToKrc(x), phi2(pp)) .* xToKrc_weights(x);
    % [~, x_weights_Ephi(:, pp)] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
    %     WeightingFunction=w);

    for cc = 0:ccOrder
        moments(cc + 1, 1) = 0.5 ./ ccL .* integral(...
            @(x) w(x) ...
            .* (-0 + cos(2*cc * acot(sqrt(x ./ ccL)))), ...
            0, inf) - 0*2.4570;
    end
    weights = ifft([moments; moments(end - 1:-1:2)]);
    weights = weights(1:ccOrder + 1);
    weights(2:end - 1) = 2*weights(2:end - 1);

    x_weights_Ephi(:, pp) = 2*ccL .* weights(2:end);
end

x_weights_E2 = sum(x_weights_Ephi .* phi2_weights, 2)


errdB2 = db(x_weights_E1 - x_weights_E2)

%% Test Spec
[rComp1] = fftCoordinates(krComp1, ApplyFftShift=true);
[rComp2] = fftCoordinates(krComp2, ApplyFftShift=true);

specTest1 = 1j * krComp1 .* modeFunE_kr(krComp1);
spatTest1 = fftshift(ifft(ifftshift(specTest1)));

specTest2 = 1j * krComp2 .* modeFunE_kr(krComp2);
spatTest2 = interp1(rComp1, spatTest1, rComp2, "spline");

specTest3 = fftshift(fft(ifftshift(spatTest2))) .* numel(rComp1) ./ numel(rComp2);

figure;
plot(rComp1, real(spatTest1), "", LineWidth=1.5);
grid on;

figure;
plot(krComp1, db(specTest1), "", LineWidth=1.5);
grid on;

hold on;
plot(krComp2, db(specTest2), "", LineWidth=1.5);
grid on;


plot(krComp2, db(specTest3), "", LineWidth=1.5);
grid on;







