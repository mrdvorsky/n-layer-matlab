function [krc, krc_weights_H] = getContourWeights_circularTE(N, Nrho, Nphi, L, Lc, wgR, krToKrc, krToKrcPrime)
%GETCONTOURWEIGHTS Helper functions for contour integration.
%
% Author: Matt Dvorsky

arguments
    N(1, 1);
    Nrho(1, 1);
    Nphi(1, 1);
    L(1, 1);
    Lc(1, 1);
    wgR(1, 1);
    krToKrc;
    krToKrcPrime;
end

%% Get Spectrum Functions
[~, specFunEy] = nLayerOpenEnded.getSpectrumCircular(...
    wgR, 0, 1, "TE");

%% Helper
modeFunH = @(kr) 2*pi * kr .* sum(specFunEy(kr, 0, kr, 0).^2, 3);

%% Compute Weights and Nodes
MkHfun = @(kr) modeFunH(krToKrc(kr)) .* krToKrcPrime(kr) .* sqrt(1 + krToKrc(kr).^2);

[x(:, 1), w(:, 1)] = fejer2(Nrho, 0, pi);
MkH = (2*L) * sum(w .* MkHfun(L * cot(0.5*x).^2) ./ (1 - cos(x)).^2 .* sin((1:N).*x), 1);

[kr, krc_weights_H(:, 1)] = fejer2_halfOpen(N, L, ...
    WeightingMoments=MkH);

krc = krToKrc(kr);

krc_weights_H = krc_weights_H ./ sqrt(1 + (krc).^2);


end

