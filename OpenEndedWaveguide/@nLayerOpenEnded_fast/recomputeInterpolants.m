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
    O.modeStructs = O.defineWaveguideModes();
end

%% Set Excitation/Receive Modes
if any([O.modeStructs.IsExcitationMode]) && any([O.modeStructs.IsReceiveMode])
    O.excitationModeIndices = find([O.modeStructs.IsExcitationMode]);
    O.receiveModeIndices = find([O.modeStructs.IsReceiveMode]);
else
    for ii = O.excitationModeIndices
        O.modeStructs(ii).IsExcitationMode = true;
    end
    for ii = O.receiveModeIndices
        O.modeStructs(ii).IsReceiveMode = true;
    end
end

%% Update Mode Types, Counts, and Cutoffs
O.modeTypes = [O.modeStructs.ModeType];
O.cutoffWavenumbers = O.modeStructs.CutoffWavenumber;

O.numModes_TE = sum(strcmp(O.modeTypes, "TE"));
O.numModes_TM = sum(strcmp(O.modeTypes, "TM"));
O.numModes_Hybrid = sum(strcmp(O.modeTypes, "Hybrid"));
O.numModes = O.numModes_TE + O.numModes_TM + O.numModes_Hybrid;

%% Fixed Point Integration Weights and Nodes
% Compute AhHat and AeHat at kRhoP coordinates.
[O.fixed_kr, O.fixed_Ah, O.fixed_Ae] = ...
    O.computeAhat();


O.fixed_Ah = reshape(O.fixed_Ah, numel(O.fixed_kr), 1, []);
O.fixed_Ae = reshape(O.fixed_Ae, numel(O.fixed_kr), 1, []);

% O.fixed_kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
% O.fixed_AhHat = AhHat .* weights;
% O.fixed_AeHat = AeHat .* weights;

end





