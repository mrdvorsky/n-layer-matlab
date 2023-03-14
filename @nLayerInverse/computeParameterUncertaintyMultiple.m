function [Uncertainty] = computeParameterUncertaintyMultiple(NLsolver, NL, f, options)
%COMPUTEPARAMETERUNCERTAINTYMULTIPLE Calculates uncertainty in parameter estimates.
% Computes the uncertainty in the structure values that are solved for
% using the "solveStructureMultiple" function. This function works similar
% to the "solveStructureMultiple" function, except it takes triplets of
% an nLayerInverse object, nLayerForward object, and a frequency vector
% (with no measurements).
%
% Example Usage (for multi-thickness, multi-band MUT measurement):
%   NLsolver1.setInitialValues(Er=[1, 4-0.1j, 1], Thk=[10, 2,  10]);
%   NLsolver2.setInitialValues(Er=[1, 4-0.1j, 1], Thk=[40, 10, 40]);
%   NLsolver1.setLayersToSolve(Erp=[2], Erpp=[2]);
%   NLsolver2.setLayersToSolve(Erp=[2], Erpp=[2]);
%   [Uncert] = nLayerInverse.computeParameterUncertaintyMultiple(...
%       NLsolver1, NL1, f1, ...
%       NLsolver2, NL2, f2, ...
%       NoiseStd=0.03);
%
% Example Usage (for multi-standoff open-ended measurements):
%   NLsolver1.setInitialValues(Er=[1, 2-0.01j], Thk=[0,  20]);
%   NLsolver2.setInitialValues(Er=[1, 2-0.01j], Thk=[10, 20]);
%   NLsolver1.setLayersToSolve(Erp=[2], Erpp=[2], Thk=[2]);
%   NLsolver2.setLayersToSolve(Erp=[2], Erpp=[2], Thk=[2]);
%   [Uncert] = nLayerInverse.computeParameterUncertaintyMultiple(...
%       NLsolver1, NL, f, ...
%       NLsolver2, NL, f, ...
%       NoiseStd=0.005);
%
% This function is mostly useful for predicting the parameter uncertainty
% given a specific multilayered structure measurement setup. The
% uncertainty is calculated assuming a specific uncertainty in the
% measurement parameters, specified using the "NoiseStd" named parameter,
% which has a default value of -40 dB or 0.01.
%
% Inputs:
%   NLsolver (Repeating) - A valid nLayerInverse object. Each must have
%       the same number of parameters to solve.
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

arguments(Repeating)
    NLsolver(1, 1) {mustBeA(NLsolver, "nLayerInverse")};
    NL(1, 1) {mustBeA(NL, "nLayerForward")};
    f(:, 1) {mustBeNonempty};
end

arguments
    options.NoiseStd(1, 1) {mustBeNonnegative} = 0.01;
end

%% Construct Linearized Ranges and Initial Guesses
[xInitial, ~, ~] = NLsolver{1}.constructInitialValuesAndRanges();

[~, gam] = nLayerInverse.calculateError(NLsolver, xInitial, NL, f, num2cell(zeros(length(NL), 1)));

%% Create Error Function
errorFunctionVector = @(x) nLayerInverse.calculateError(NLsolver, x, NL, f, gam);

%% Calculate Jacobian
[~, ~, ~, ~, ~, ~, J]  = lsqnonlin(errorFunctionVector, xInitial, [], [], ...
    optimoptions("lsqnonlin", Display="none"));

%% Compute Uncertainty from Jacobian
hess = pinv(full(J).' * full(J), 0) .* options.NoiseStd.^2;
xUncertainty = sqrt(diag(hess));

Uncertainty = cell(numel(NLsolver), 1);
for ii = 1:numel(Uncertainty)
    [er, ur, thk] = NLsolver{ii}.extractStructure(xInitial, f);
    [Uncertainty{ii}.erLower, Uncertainty{ii}.urLower, Uncertainty{ii}.thkLower] = ...
        NLsolver{ii}.extractStructure(xInitial - xUncertainty, f);
    [Uncertainty{ii}.erUpper, Uncertainty{ii}.urUpper, Uncertainty{ii}.thkUpper] = ...
        NLsolver{ii}.extractStructure(xInitial + xUncertainty, f);

    Uncertainty{ii}.er  = complex(...
        max(abs(real(er - Uncertainty{ii}.erLower)), abs(real(er - Uncertainty{ii}.erUpper))), ...
        max(abs(imag(er - Uncertainty{ii}.erLower)), abs(imag(er - Uncertainty{ii}.erUpper))));
    Uncertainty{ii}.ur  = complex(...
        max(abs(real(ur - Uncertainty{ii}.urLower)), abs(real(ur - Uncertainty{ii}.urUpper))), ...
        max(abs(imag(ur - Uncertainty{ii}.urLower)), abs(imag(ur - Uncertainty{ii}.urUpper))));
    Uncertainty{ii}.thk = ...
        max(abs(thk - Uncertainty{ii}.thkLower), abs(thk - Uncertainty{ii}.thkUpper));
end

end

