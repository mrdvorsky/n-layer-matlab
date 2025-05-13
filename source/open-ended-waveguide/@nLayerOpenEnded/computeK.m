function [K] = computeK(self, f)
%Computes the matrix K.
% This function computes the matrix K as a function of each frequency
% specified by "f", which is used to compute the mode S-parameter matrix.
%
% Example Usage:
%   [A] = self.computeA(f, er, ur, thk);
%   [K] = self.computeK(f);
%   idMat = eye(size(A, 1));
%   Smn = pagemldivide(idMat + A.*K, idMat - A.*K);
%
% Author: Matt Dvorsky

arguments
    self nLayerOpenEnded;

    f(1, 1, :) double;
end

%% Check Frequency Range
if min(f) < self.frequencyRange(1) || max(f) > self.frequencyRange(2)
    error("nLayerOpenEnded:outsideOfFrequencyRange", ...
        "At least one frequency that was passed in is outside " + ...
        "the operating range specified [%g, %g], and so the results " + ...
        "may not be accurate. Please modify the 'frequencyRange' " + ...
        "parameter.", self.frequencyRange(1), self.frequencyRange(2));
end

%% Get Waveguide Fill er and ur for Each Frequency
wgEr = self.waveguideEr;
if ~isnumeric(wgEr)     % If wgEr is a function Handle.
    wgEr = wgEr(f);
end

wgUr = self.waveguideEr;
if ~isnumeric(wgUr)     % If wgUr is a function Handle.
    wgUr = wgUr(f);
end

%% Mode Coefficients
k0 = 2*pi .* f ./ self.speedOfLight;

beta = sqrt(k0.^2 .* wgEr .* wgUr - self.mode_kc0(:).^2);
beta = complex(real(beta), -abs(imag(beta)));

%% Compute K
isTE = strcmp(self.modeTypes, "TE");

beta_TE_over_k0 = (wgUr .* k0) ./ beta;
beta_TM_over_k0 = beta ./ (wgEr .* k0);
betaAll = beta_TM_over_k0;
betaAll(isTE, :) = beta_TE_over_k0(isTE, :);

K = sqrt(betaAll);
K = (K) .* (pagetranspose(K));

end

