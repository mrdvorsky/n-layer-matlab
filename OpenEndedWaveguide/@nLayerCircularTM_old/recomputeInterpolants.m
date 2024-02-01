function [] = recomputeInterpolants(O)
%RECOMPUTEINTERPOLANTS Recompute interpolation functions, structs, etc.
% This function should be called anytime any critical parameter is changed,
% and must be called before calling calculate. This function is
% automatically called after creating an nLayerCircularTM object.
%
% List of critical parameters:
%   waveguideR;
%   speedOfLight;
%   modesTM;
%   interpolationPoints_kRho;
%   integralPointsFixed_kRho;
%   integralInitialSegmentCount;
%
% Example Usage:
%   NL = nLayerRectangular(...);
%   NL.*criticalParameter* = *newValue*;
%   NL.recomputeInterpolants();
%
% Author: Matt Dvorsky

%% Calculate Mode Cutoffs
O.numModes = numel(O.modesTM);

O.modeCutoffs = zeros(O.numModes, 1);
for ii = 1:O.numModes
    % Use zeros of 0th order bessel function to calculate mode cutoff
    % wavenumbers. A reasonable approximation for the nth zero is 
    % pi*(n - 0.25). This is used as the initial guess.
    O.modeCutoffs(ii) = fzero(@(x) besselj(0, x), ...
        pi*(O.modesTM(ii) - 0.25)) ./ O.waveguideR;
end

% Scale factor for change of variables between tau and tauP
O.integralScaleFactor = (2*pi) ./ O.waveguideR;

%% Compute the Matrix A at Various Values of kRhoP
% Compute Ah and interpolation lookup tables as a function of kRhoP.
kRhoP(:, 1) = linspace(0, 1, O.interpolationPoints_kRho);

% Compute AeHat at kRhoP coordinates.
[AeHat] = O.computeAhat(kRhoP);

% Store in a table that is used in the "integrandAhat" function.
O.table_AeHat = AeHat;

%% Fixed Point Integration Weights and Nodes
% For lossy structures, generally no adaptive meshing is needed. In those
% cases we can use a precomputed set of weights and nodes, instead of
% computing them on the fly. This is generally 3 to 4 times faster than
% when using the adaptive integration.
[kRhoP, weights, errWeights] = fejer2(O.integralPointsFixed_kRho, 0, 1);

% Compute AeHat at kRhoP coordinates.
[AeHat] = O.computeAhat(kRhoP);

% Store computed matrices. Also, precompute kRho using kRhoP. These are
% used in the "computeA" function.
O.fixed_kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
O.fixed_AeHat = AeHat .* weights;
O.fixed_errorAeHat = AeHat .* errWeights;

%% First Integration Pass Precomputation
% The initial pass of the adaptive integral algorithm always uses the same
% nodes (i.e., the same evaluation coordinates of kRho). Thus, we can 
% precompute the values of the integrand for A at those coordinates,
% instead of needing to perform an interpolation. These are used in the
% "integrandAhat" function.
[kRhoP, ~, ~] = nLayer.gaussKronrod(...
    O.integralInitialSegmentCount, 0, 1);

% Compute AeHat at kRhoP coordinates.
[AeHat] = O.computeAhat(kRhoP);

% Store computed matrices. Also, precompute kRho using kRhoP. These are
% used in the "integrandAhat" function.
O.init_kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
O.init_AeHat = AeHat;

end

