function [krc, krc_weights_E, krc_weights_H, krc_error_E, krc_error_H] = getContourWeights_customPath(N, Nrho, Nphi, L, wgA, wgB, krToKrc, krToKrcPrime)
%GETCONTOURWEIGHTS Helper functions for contour integration.
%
% Author: Matt Dvorsky

arguments
    N(1, 1);
    Nrho(1, 1);
    Nphi(1, 1);
    L(1, 1);
    wgA(1, 1);
    wgB(1, 1);
    krToKrc;
    krToKrcPrime;
end

%% Get Spectrum Functions
% [~, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
%     wgA, wgB, 1, 0, "TE");
NL = nLayerRectangular(1, 0, waveguideA=wgA, waveguideB=wgB);
specFunEy = NL.modeStructs.EySpec;

%% Helper
[phi(1, 1, :), phi_weights(1, 1, :)] = fejer2(Nphi, 0, pi/2);
phi_weights = 4*phi_weights;

modeFunE = @(kr) kr .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* sin(phi).^2 .* phi_weights, 3);
modeFunH = @(kr) kr .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* cos(phi).^2 .* phi_weights, 3);

%% Compute Weights and Nodes
MkEfun = @(kr) modeFunE(krToKrc(kr)) .* krToKrcPrime(kr) ./ sqrt(1 + krToKrc(kr).^2);
MkHfun = @(kr) modeFunH(krToKrc(kr)) .* krToKrcPrime(kr) .* sqrt(1 + krToKrc(kr).^2);

% MkEfun = @(kr) modeFunE(krToKrc(kr)) .* krToKrcPrime(kr);
% MkHfun = @(kr) modeFunH(krToKrc(kr)) .* krToKrcPrime(kr);

[x(:, 1), w(:, 1)] = trap(Nrho, 0, pi);
MkH = (2*L) * sum(w .* MkHfun(L * cot(0.5*x).^2) ./ (1 - cos(x)).^2 .* sin((1:N).*x), 1);
MkE = (2*L) * sum(w .* MkEfun(L * cot(0.5*x).^2) ./ (1 - cos(x)).^2 .* sin((1:N).*x), 1);

[kr, krc_weights_E(:, 1), krc_error_E(:, 1)] = fejer2_halfOpen(N, L, ...
    WeightingMoments=MkE);

[~, krc_weights_H(:, 1), krc_error_H(:, 1)] = fejer2_halfOpen(N, L, ...
    WeightingMoments=MkH);

krc = krToKrc(kr);

krc_weights_E = krc_weights_E .* sqrt(1 + (krc).^2);
krc_weights_H = krc_weights_H ./ sqrt(1 + (krc).^2);

krc_error_E = krc_error_E .* sqrt(1 + (krc).^2);
krc_error_H = krc_error_H ./ sqrt(1 + (krc).^2);


end

