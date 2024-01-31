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
if isempty(O.modeStructs)
    modeStruct = O.defineWaveguideModes();
else
    modeStruct = O.modeStructs{1};
end

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

%% Fixed Point Integration Weights and Nodes
% Compute AhHat and AeHat at kRhoP coordinates.
[O.fixed_kRho, O.fixed_Ah_weights, O.fixed_Ae_weights] = ...
    O.computeAhat(modeStruct);


O.fixed_Ah_weights = reshape(O.fixed_Ah_weights, numel(O.fixed_kRho), 1, []);
O.fixed_Ae_weights = reshape(O.fixed_Ae_weights, numel(O.fixed_kRho), 1, []);

% O.fixed_kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
% O.fixed_AhHat = AhHat .* weights;
% O.fixed_AeHat = AeHat .* weights;

end





