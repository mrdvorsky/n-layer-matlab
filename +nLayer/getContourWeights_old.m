function [krc, krc_weights_E, krc_weights_H] = nlGetContourWeights(N, Nrho, Nphi, L, Lc, wgA, wgB)
%GETCONTOURWEIGHTS Helper functions for contour integration.
%
% Author: Matt Dvorsky

arguments
    N(1, 1);
    Nrho(1, 1);
    Nphi(1, 1);
    L(1, 1);
    Lc(1, 1);
    wgA(1, 1);
    wgB(1, 1);
end

%% Get Spectrum Functions
[~, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    wgA, wgB, 1, 0, "TE");

%% Helper
[phi(1, 1, :), phi_weights(1, 1, :)] = fejer2(Nphi, 0, pi/2);
phi_weights = 4*phi_weights;

modeFunE = @(kr) kr .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* sin(phi).^2 .* phi_weights, 3);
modeFunH = @(kr) kr .* sum(specFunEy(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* cos(phi).^2 .* phi_weights, 3);

%% Compute Weights and Nodes
[krc, momentKernel] = getContourIntegrals(N, L, Lc);
[kr(:, 1), kr_weights(:, 1)] = fejer2_halfOpen(Nrho, L);

MkE = sum((kr_weights ./ sqrt(1 + kr.^2) .* modeFunE(kr) .* momentKernel(kr, 1:N)), 1);
MkH = sum((kr_weights .* sqrt(1 + kr.^2) .* modeFunH(kr) .* momentKernel(kr, 1:N)), 1);

[~, krc_weights_E(:, 1)] = fejer2_halfOpen(N, L, ...
    WeightingMoments=MkE);
krc_weights_E = krc_weights_E .* sqrt(1 + (krc).^2);

[~, krc_weights_H(:, 1)] = fejer2_halfOpen(N, L, ...
    WeightingMoments=MkH);
krc_weights_H = krc_weights_H ./ sqrt(1 + (krc).^2);


end

