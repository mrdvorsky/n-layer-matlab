function [AeHat] = computeAhat(O, kRhoP)
%COMPUTEAHAT Computes the matrix AeHat(kRho).
%   This function computes the matrix AeHat(kRho) as a function of kRho.
%   The outputs of this function can be used to compute the matrix Ae by
%   integrating over kRhoP over the interval [0, 1]. Note that the change
%   of variables from kRhoP to kRho will be applied in this function
%   automatically.
%
% Example Usage:
%   [AeHat] = O.computeAhat(kRhoP);
%
% Inputs:
%   kRhoP - A vector of kRhoP coordinates (coordinate transform of kRho).
%       Values should be in the interval [0, 1].
% Outputs:
%   AeHat - Matrix computed as a function of kRhoP that can be used to
%       compute the matrix A. The size will be numel(kRhoP) by O.numModes
%       by O.numModes.
%
% Author: Matt Dvorsky

arguments
    O;
    kRhoP(:, 1);
end

%% Integral Change of Variables
% The integral needs to be evaluated from kRho = [0, inf). However, a change
% of variables kRho = L(1 - kRhoP)/kRhoP is used here so that the interpolant
% can be uniform in (0, 1].
L = O.integralScaleFactor;
kRho = L * (1 - kRhoP) ./ kRhoP;

% Weighting function to account for change of variables.
weights_kRho = L ./ (kRhoP.^2);

%% Compute Waveguide Mode Cutoffs
kc0i(1, 1, :, 1) = O.modeCutoffs;
kc0j(1, :, 1, 1) = O.modeCutoffs;

%% Compute AeHat(kRhoP)
AeHat = weights_kRho .* kc0j .* besselj(1, O.waveguideR .* kc0j) ...
    .* kRho.^3 .* besselj(0, O.waveguideR * kRho).^2 ...
    ./ ((kRho.^2 - kc0i.^2) .* (kRho.^2 - kc0j.^2));

%% Fix Nans Caused by Singularities At Endpoints
AeHat(kRhoP == 0, :, :) = 0;
AeHat(kRhoP == 1, :, :) = 0;

end

