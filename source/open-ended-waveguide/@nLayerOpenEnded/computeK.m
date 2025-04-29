function [K] = computeK(O, f)
%COMPUTEK Computes the matrix K.
% This function computes the matrix K as a function of each frequency
% specified by "f", which is used to compute the mode S-parameter matrix.
%
% Example Usage:
%   [A] = O.computeA(f, er, ur, thk);
%   [K] = O.computeK(f);
%   idMat = eye(size(A, 1));
%   Smn = pagemldivide(idMat + A.*K, idMat - A.*K);
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
wgEr = O.waveguideEr;
if ~isnumeric(wgEr)     % If wgEr is a function Handle.
    wgEr = wgEr(f);
end

wgUr = O.waveguideEr;
if ~isnumeric(wgUr)     % If wgUr is a function Handle.
    wgUr = wgUr(f);
end

%% Mode Coefficients
k0 = 2*pi .* f ./ O.speedOfLight;

beta = sqrt(k0.^2 .* wgEr .* wgUr - O.mode_kc0(:).^2);
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

