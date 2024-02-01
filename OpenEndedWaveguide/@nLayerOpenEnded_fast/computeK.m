function [K] = computeK(O, f)
%COMPUTEK Computes kA and kB.
% Computes the matrices kA and kB, which are used to compute the
% unnormalized mode S-parameter matrix.
%
% Inputs:
%   f - vector of frequencies to consider.
% Outputs:
%   kA, kB - Matrices used to compute S-parameters.
%
% Example Usage (for single frequency, unnormalized S-parameter matrix):
%   [A] = O.computeA(f(ii), er, ur, thk);
%   [B] = O.computeB();
%   [kA, kB] = O.computeK(f(ii));
%   S = (A.*kA + B.*kB) \ (-A.*kA + B.*kB);
%
% Although the example above shows usage with a scalar value for "f", the
% input "f" can be a vector. In this case, the size of the 3rd dimension
% of each output matrix will be equal to numel(f).
%
% Author: Matt Dvorsky

arguments
    O;
    f(1, 1, :) double;
end

%% Mode Coefficients
k0 = 2*pi .* f ./ O.speedOfLight;

wgEr = cell2mat(cellfun(@(fun) fun(k0), O.waveguideEr, UniformOutput=false));
wgUr = cell2mat(cellfun(@(fun) fun(k0), O.waveguideUr, UniformOutput=false));

beta = conj(sqrt(k0.^2 .* wgEr .* wgUr - O.cutoffWavenumbers.^2));
beta = complex(real(beta), -abs(imag(beta)));

%% Compute K
isTE = strcmp(O.modeTypes, "TE");

beta_TE_over_k0 = (wgUr .* k0) ./ beta;
beta_TM_over_k0 = beta ./ (wgEr .* k0);
betaAll = beta_TM_over_k0;
betaAll(isTE, :) = beta_TE_over_k0(isTE, :);

K = sqrt(betaAll);
K = (K) .* (pagetranspose(K));

end

