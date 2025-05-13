function [Uncertainty] = computeParameterUncertainty(self, NL, f, options)
%Calculates uncertainty in parameter estimates.
% Computes the uncertainty in the structure values that are solved for
% using the "solveStructure" function. This function works similar to the
% "solveStructure" function, except it takes pairs of an nLayerForward
% object and a frequency vector (with no measurements).
%
% Example Usage (simple case, 2 layer infinite halfspace):
%   NLsolver.setInitialValues(Er=[1, 4-1j], Thk=[10, inf]);
%   NLsolver.setLayersToSolve(Er=[2]);
%   [Uncert] = NLsolver.computeParameterUncertainty(NL, f, gam, ...
%       NoiseStd=0.03);
%
% Example Usage (for multi-band open-ended measurements):
%   NLsolver.setInitialValues(Er=[1, 4-1j], Thk=[1, 10]);
%   NLsolver.setLayersToSolve(Erp=[2], Erpp=[2], Thk=[1, 2]);
%   [Uncert] = NLsolver.computeParameterUncertainty(...
%       NL1, f1, gam1, ...
%       NL2, f2, gam2, ...
%       NoiseStd=0.03);
%
%
% This function is mostly useful for predicting the parameter uncertainty
% given a specific multilayered structure measurement setup. The
% uncertainty is calculated assuming a specific uncertainty in the
% measurement parameters, specified using the "NoiseStd" named parameter,
% which has a default value of -40 dB or 0.01.
%
% Inputs:
%   NL (Repeating) - A valid nLayerForward object.
%   f (Repeating) - Vector of frequencies to pass to NL.
%
% Outputs:
%   Uncertainty - Cell array of structs containing the calculated output
%       parmeter uncertainties for each input set.
%
% Named Arguments:
%   NoiseStd (0.01) - Uncertainty value to use for the measurement data.
%
% Author: Matt Dvorsky

arguments
    self nLayerInverse;
end
arguments(Repeating)
    NL(1, 1) nLayerForward;
    f(:, 1) {mustBeNonempty};
end
arguments
    options.NoiseStd(1, 1) {mustBeNonnegative} = 0.01;
end

%% Calculate Uncertainty Using "computeParameterUncertaintyMultiple"
inputParams = [repmat({self}, 1, numel(NL)); NL; f];
[Uncertainty] = nLayerInverse.computeParameterUncertaintyMultiple(...
    inputParams{:}, NoiseStd=options.NoiseStd);

Uncertainty = Uncertainty{1};   % Get first element of Uncertainty.

end

