function [] = recomputeInterpolants(O)
%RECOMPUTEINTERPOLANTS Recompute interpolation functions, structs, etc.
% This function should be called anytime any critical parameter is changed,
% and must be called before calling calculate. This function is
% automatically called after creating an nLayerRectangular object.
%
% List of critical parameters:
%   waveguideA;
%   waveguideB;
%   speedOfLight;
%   modesTE;
%   interpolationPoints_kRho;
%   integralPointsFixed_kRho;
%   integralPoints_kPhi;
%   integralInitialSegmentCount;
%
% Example Usage:
%   NL = nLayerRectangular(...);
%   NL.*criticalParameter* = *newValue*;
%   NL.recomputeInterpolants();
%
% Author: Matt Dvorsky

%% Update TM Modes
O.numModes = numel(O.modeSpectrumTE_Hx) + numel(O.modeSpectrumTM_Ex);

% Scale factor for change of variables between kRho and kRhoP
% O.integralScaleFactor = ;

%% Compute the Matrix A at Various Values of kRhoP
% Compute AhHat and AeHat interpolation lookup tables as a function of kRhoP.
kRhoP(:, 1) = linspace(0, 1, O.interpolationPoints_kRho);

% Compute AhHat and AeHat at kRhoP coordinates.
[AhHat, AeHat] = O.computeAhat(kRhoP);

% Combine AeHat and AhHat into one array for faster interpolation.
% Store in a table that is used in the "integrandAhat" function.
O.table_AheHat = cat(4, AhHat, AeHat);

%% Fixed Point Integration Weights and Nodes
% For lossy structures, generally no adaptive meshing is needed. In those
% cases we can use a precomputed set of weights and nodes, instead of
% computing them on the fly. This is generally 3 to 4 times faster than
% when using the adaptive integration.
[kRhoP, weights, errWeights] = O.fejer2(O.integralPointsFixed_kRho, 0, 1);

% Compute AhHat and AeHat at kRhoP coordinates.
[AhHat, AeHat] = O.computeAhat(kRhoP);

% Store computed matrices. Also, precompute kRho using kRhoP. These are
% used in the "computeA" function.
O.fixed_kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
O.fixed_AhHat = AhHat .* weights;
O.fixed_AeHat = AeHat .* weights;
O.fixed_errorAhHat = AhHat .* errWeights;
O.fixed_errorAeHat = AeHat .* errWeights;

%% First Integration Pass Precomputation
% The initial pass of the adaptive integral algorithm always uses the same
% nodes (i.e., the same evaluation coordinates of kRho). Thus, we can 
% precompute the values of the integrand for A at those coordinates,
% instead of needing to perform an interpolation. These are used in the
% "integrandAhat" function.
[kRhoP, ~, ~] = O.gaussKronrod(...
    O.integralInitialSegmentCount, 0, 1);

% Compute AhHat and AeHat at kRhoP coordinates.
[AhHat, AeHat] = O.computeAhat(kRhoP);

% Store computed matrices. Also, precompute kRho using kRhoP. These are
% used in the "integrandAhat" function.
O.init_kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
O.init_AhHat = AhHat;
O.init_AeHat = AeHat;

end

