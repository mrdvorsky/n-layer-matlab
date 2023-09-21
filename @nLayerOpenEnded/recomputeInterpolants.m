function [] = recomputeInterpolants(O)
%RECOMPUTEINTERPOLANTS Recompute interpolation functions, structs, etc.
% This function is called whenever a parameter changes that would change
% the mode spectrums.
%
% Author: Matt Dvorsky

arguments
    O;
end

%% Get Waveguide Mode Specifications
modeStruct = O.defineWaveguideModes();

O.outputModes_TE = modeStruct.OutputModes_TE;
O.outputModes_TM = modeStruct.OutputModes_TM;
O.outputModes_Hybrid = modeStruct.OutputModes_Hybrid;

%% Update Mode Counts and Cutoffs
O.numModes_TE = numel(modeStruct.SpecEx_TE);
O.numModes_TM = numel(modeStruct.SpecEx_TM);
O.numModes_Hybrid = numel(modeStruct.SpecEx_Hybrid);
O.numModes = O.numModes_TE + O.numModes_TM + O.numModes_Hybrid;

O.cutoffBeta_TE = modeStruct.CutoffBeta_TE;
O.cutoffBeta_TM = modeStruct.CutoffBeta_TM;
O.cutoffBeta_Hybrid = modeStruct.CutoffBeta_Hybrid;

%% Update Integral Scale Factor
O.integralScaleFactor = modeStruct.IntegralScaleFactor;

%% Compute Ahat at all possible values of kRhoP
kRhoP(:, 1) = linspace(0, 1, O.interpolationPoints_kRho);

% Compute AhHat and AeHat at kRhoP coordinates.
[AhHat, AeHat] = O.computeAhat(kRhoP, modeStruct);
modeStruct.CheckModeScalingAndOrthogonality = false;

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
[AhHat, AeHat] = O.computeAhat(kRhoP, modeStruct);

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
[AhHat, AeHat] = O.computeAhat(kRhoP, modeStruct);

% Store computed matrices. Also, precompute kRho using kRhoP. These are
% used in the "integrandAhat" function.
O.init_kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
O.init_AhHat = AhHat;
O.init_AeHat = AeHat;

end





