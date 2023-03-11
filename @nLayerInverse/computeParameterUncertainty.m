function [Uncertainty] = computeParameterUncertainty(O, NL, f, options)
%COMPUTEPARAMETERUNCERTAINTY Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
end

arguments(Repeating)
    NL;
    f(:, 1);
end

arguments
    options.NoiseStd = 0.001;
end

%% Construct Linearized Ranges and Initial Guesses
O.validate();
[xInitial, ~, ~] = O.constructInitialValuesAndRanges();

[~, gam] = O.calculateError({O}, xInitial, NL, f, num2cell(zeros(length(NL), 1)));

%% Create Error Function
errorFunctionVector = @(x) O.calculateError({O}, x, NL, f, gam);

%% Calculate Jacobian
[~, ~, ~, ~, ~, ~, J]  = lsqnonlin(errorFunctionVector, xInitial, [], [], ...
    optimoptions(O.localOptimizerOptions, Display="none"));

%% Compute Uncertainty from Jacobian
hess = inv(full(J).' * full(J)) .* 0.5 * options.NoiseStd.^2;
xUncertainty = sqrt(diag(hess));

[er, ur, thk] = O.extractStructure(xInitial, f);
[Uncertainty.erLower, Uncertainty.urLower, Uncertainty.thkLower] = ...
    O.extractStructure(xInitial - xUncertainty, f);
[Uncertainty.erUpper, Uncertainty.urUpper, Uncertainty.thkUpper] = ...
    O.extractStructure(xInitial + xUncertainty, f);

Uncertainty.er  = complex(...
    max(abs(real(er  - Uncertainty.erLower)),  abs(real(er  - Uncertainty.erUpper))), ...
    max(abs(imag(er  - Uncertainty.erLower)),  abs(imag(er  - Uncertainty.erUpper))));
Uncertainty.ur  = complex(...
    max(abs(real(ur  - Uncertainty.urLower)),  abs(real(ur  - Uncertainty.urUpper))), ...
    max(abs(imag(ur  - Uncertainty.urLower)),  abs(imag(ur  - Uncertainty.urUpper))));
Uncertainty.thk = ...
    max(abs(thk - Uncertainty.thkLower), abs(thk - Uncertainty.thkUpper));

end

