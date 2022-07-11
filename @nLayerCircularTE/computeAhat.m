function [AhHat] = computeAhat(O, kRhoP)
%COMPUTEAHAT Computes the matrix AhHat(kRho).
%   This function computes the matrix AhHat(kRho) as a function of kRho.
%   The outputs of this function can be used to compute the matrix Ah by
%   integrating over kRhoP over the interval [0, 1]. Note that the change
%   of variables from kRhoP to kRho will be applied in this function
%   automatically.
%
% Example Usage:
%   [AhHat] = O.computeAhat(kRhoP);
%
% Inputs:
%   kRhoP - A vector of kRhoP coordinates (coordinate transform of kRho).
%       Values should be in the interval [0, 1].
% Outputs:
%   AhHat - Matrix computed as a function of kRhoP that can be used to
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

%% Compute AhHat(kRhoP)
AhHat = weights_kRho .* kc0j.^2 .* kc0i .* besselj(0, O.waveguideR .* kc0j) ...
    .* kRho .* besselj(1, O.waveguideR * kRho).^2 ...
    ./ ((kRho.^2 - kc0i.^2) .* (kRho.^2 - kc0j.^2));

%% Fix Nans Caused by Singularities At Endpoints
AhHat(kRhoP == 0, :, :) = 0;
AhHat(kRhoP == 1, :, :) = 0;

end

