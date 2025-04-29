clc;
clear;
close all;

%% Inputs
Lc = 1;
ccL = 1;

ccOrder = 10;

N1 = 400;
N2 = 40;

krMax = 100;

%% nLayer
NL = nLayerRectangular(1, 0, waveguideBand="ka", convergenceAbsTol=1e-7);

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

%% Helper
phi1(1, 1, 1, :) = 2*pi * ((1:4*N1) - 0.5) ./ (4*N1);
phi1_weights = 0*phi1(1:N1) + 4*(2*pi ./ (4*N1));
phi1 = phi1(1:N1);

modeFunE_kr = @(kr) sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
    .* sin(phi1).^2 .* phi1_weights, 4);

alpha = @(x) x + 1j .* Lc .* x ./ (Lc + x);
alphaP = @(x) 1 + 1j .* Lc.^2 ./ (x + Lc).^2;

beta = @(x) 0.5 * (sqrt(x.^2 + (2-2j)*x + 2j) + x - 1 - 1j);
betaP = @(x) 0.5 * (1 + (x + 1 - 1j) ...
    ./ sqrt(x.^2 + (2-2j)*x + 2j));

% beta = @(x) x + 1j .* Lc .* x ./ (Lc + x);
% betaP = @(x) 1 + 1j .* Lc.^2 ./ (x + Lc).^2;
% 
% alpha = @(x) 0.5 * (sqrt(x.^2 + (2-2j)*x + 2j) + x - 1 - 1j);
% alphaP = @(x) 0.5 * (1 + (x + 1 - 1j) ...
%     ./ sqrt(x.^2 + (2-2j)*x + 2j));

%% Helper 2
% phi2(1, :) = 2*pi * ((1:4*N2) - 0.5) ./ (4*N2);
% phi2_weights = 0*phi2 + 1*(2*pi ./ (4*N2));
% % phi2 = phi2(1:N2);

[phi2(1, :), phi2_weights(1, :)] = fejer2(N2, 0, pi/2);
phi2_weights = 4*phi2_weights;

modeFunE_krphi = @(kr, ph) specFunEy(kr .* cos(ph), kr .* sin(ph)).^2 ...
    .* sin(ph).^2;

%% Compute Weighting 1
tic;
[x, x_weights_E1] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
    WeightingFunction=@(x) modeFunE_kr(alpha(x)) .* alphaP(x));
toc;

%% Compute Weighting 2
tic;
for pp = 1:numel(phi2)
    w = @(x) modeFunE_krphi(alpha(x), phi2(pp)) .* alphaP(x);

    for cc = 0:ccOrder
        moments(cc + 1, 1) = 0.5 ./ ccL .* integral(...
            @(x) w(x) ...
            .* (cos(2*cc * acot(sqrt(x ./ ccL)))), ...
            0, inf);
    end
    weights = ifft([moments; moments(end - 1:-1:2)]);
    weights = weights(1:ccOrder + 1);
    weights(2:end - 1) = 2*weights(2:end - 1);

    x_weights_Ephi(:, pp) = 2*ccL .* weights(2:end);
end
toc;

x_weights_E2 = sum(x_weights_Ephi .* phi2_weights, 2)


errdB2 = db(x_weights_E1 - x_weights_E2)


%% Compute Weighting 3
tic;
for pp = 1:numel(phi2)
    w = @(x) modeFunE_krphi(alpha(x), phi2(pp)) .* alphaP(x);

    for cc = 0:ccOrder
        w2 = @(x) w(x) .* (cos(2*cc * acot(sqrt(x ./ ccL))));
        moments(cc + 1, 1) = 0.5 ./ ccL .* integral(...
            @(x) w2(beta(x)) .* betaP(x), 0, inf);
    end
    weights = ifft([moments; moments(end - 1:-1:2)]);
    weights = weights(1:ccOrder + 1);
    weights(2:end - 1) = 2*weights(2:end - 1);

    x_weights_Ephi(:, pp) = 2*ccL .* weights(2:end);
end
toc;

x_weights_E3 = sum(x_weights_Ephi .* phi2_weights, 2)


errdB3 = db(x_weights_E1 - x_weights_E3)





