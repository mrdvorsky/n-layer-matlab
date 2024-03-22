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

%% Check Frequency Range
if min(f) < O.frequencyRange(1) || max(f) > O.frequencyRange(2)
    error("At least one frequency that was passed in is outside " + ...
        "the operating range specified [%g, %g], and so the results " + ...
        "may not be accurate. Please modify the 'frequencyRange' " + ...
        "parameter.", O.frequencyRange(1), O.frequencyRange(2));
end

%% Get Waveguide Fill er and ur for Each Frequency
wgEr = zeros(O.numModes, 1, numel(f));
for ii = 1:size(wgEr, 1)
    if isnumeric(O.modeStructs(ii).WaveguideEr)
        wgEr(ii, 1, :) = O.modeStructs(ii).WaveguideEr;
    else
        wgEr(ii, 1, :) = O.modeStructs(ii).WaveguideEr(f);
    end
end

wgUr = zeros(size(wgEr));
for ii = 1:size(wgUr, 1)
    if isnumeric(O.modeStructs(ii).WaveguideUr)
        wgUr(ii, 1, :) = O.modeStructs(ii).WaveguideUr;
    else
        wgUr(ii, 1, :) = O.modeStructs(ii).WaveguideUr(f);
    end
end

%% Mode Coefficients
k0 = 2*pi .* f ./ O.speedOfLight;

beta = conj(sqrt(k0.^2 .* wgEr .* wgUr - O.mode_kc0(:).^2));
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

