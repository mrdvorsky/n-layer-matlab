function [Uncertainty] = computeParameterUncertaintyMultiple(NLsolver, NL, f, options)
%Calculates uncertainty in parameter estimates.
% Computes the uncertainty in the structure values that are solved for
% using the "solveStructureMultiple" function. This function works similar
% to the "solveStructureMultiple" function, except it takes triplets of
% an nLayerInverse object, nLayerForward object, and a frequency vector
% (with no measurements).
%
% Example Usage (for multi-thickness, multi-band MUT measurement):
%   NLsolver1.setInitialValues(Er=[1, 4-0.1j, 1], Thk=[10, 2,  10]);
%   NLsolver2.setInitialValues(Er=[1, 4-0.1j, 1], Thk=[40, 10, 40]);
%   NLsolver1.setLayersToSolve(Er=[2]);
%   NLsolver2.setLayersToSolve(Er=[2]);
%   [Uncert] = nLayerInverse.computeParameterUncertaintyMultiple(...
%       NLsolver1, NL1, f1, ...
%       NLsolver2, NL2, f2, ...
%       NoiseStd=0.03);
%
% Example Usage (for multi-standoff open-ended measurements):
%   NLsolver1.setInitialValues(Er=[1, 2-0.01j], Thk=[0,  20]);
%   NLsolver2.setInitialValues(Er=[1, 2-0.01j], Thk=[10, 20]);
%   NLsolver1.setLayersToSolve(Er=[2], Thk=[2]);
%   NLsolver2.setLayersToSolve(Er=[2], Thk=[2]);
%   [Uncert] = nLayerInverse.computeParameterUncertaintyMultiple(...
%       NLsolver1, NL, f, ...
%       NLsolver2, NL, f, ...
%       NoiseStd=0.005);
%
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
    NLsolver(1, 1) nLayerInverse;
    NL(1, 1) nLayerForward;
    f(:, 1) {mustBeNonempty};
end
arguments
    options.NoiseStd(1, 1) {mustBeNonnegative} = 0.01;
end

%% Construct Linearized Ranges and Initial Guesses
[xInitial, ~, ~, ~, ~, xAeq, xbeq] = NLsolver{1}.constructInitialValuesAndRanges();
[~, gam] = nLayerInverse.calculateError(NLsolver, xInitial, NL, f, num2cell(zeros(length(NL), 1)));

%% Create Error Function
errorFunctionVector = @(x) nLayerInverse.calculateError(NLsolver, x, NL, f, gam);

%% Calculate Jacobian
J = 0;
if ~isempty(xInitial)
    [~, ~, ~, ~, ~, ~, J]  = lsqnonlin(errorFunctionVector, xInitial, [], [], ...
        [], [], [], [], [], optimoptions("lsqnonlin", Display="iter"));
end

%% Compute Uncertainty from Jacobian
if numel(xAeq) == 0
    xUncertainty = sqrt(diag(pinv(full(J).' * full(J), 0))) .* options.NoiseStd;
else    % Take equality constraints into account.
    [R, xSubInd(:, 1)] = rref([xAeq, xbeq]);
    yInd(:, 1) = setdiff((1:numel(xInitial)).', xSubInd);

    % Write x in the form of C*y + c0.
    C = zeros(numel(xInitial), numel(yInd));
    C(yInd, :) = eye(numel(yInd));
    C(xSubInd, :) = -R(:, yInd);

    % Jacobian for reduced order equation system.
    yUncertainty = sqrt(diag(pinv(...
        C.' * (full(J).' * full(J)) * C, 0))) .* options.NoiseStd;
    xUncertainty = abs(C*yUncertainty);
end

Uncertainty = cell(numel(NLsolver), 1);
for ii = 1:numel(Uncertainty)
    [er, ur, thk] = NLsolver{ii}.extractStructure(xInitial, f);
    Uncertainty{ii}.J = J;
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

