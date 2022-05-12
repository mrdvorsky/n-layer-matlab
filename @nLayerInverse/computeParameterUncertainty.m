function [xUncertainty] = computeParameterUncertainty(O, NL, f, options)
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
[xInitial, ~, ~] = O.constructInitialValuesAndRanges();

[~, gam] = O.calculateError(xInitial, NL, f, num2cell(zeros(length(NL), 1)));

%% Create Error Function
errorFunctionVector = @(x) O.calculateError(x, NL, f, gam);

%% Calculate Jacobian
[~, ~, ~, ~, ~, ~, J]  = lsqnonlin(errorFunctionVector, xInitial, [], [], ...
    optimoptions(O.localOptimizerOptions, Display="none"));

% svd(full(J).' * full(J))
xUncertainty = inv(full(J).' * full(J)) .* 0.5 * options.NoiseStd.^2;
% xUncertainty = (pinv(full(J)) * pinv(full(J)).') .* 0.5 .* options.NoiseStd.^2;

disp(sqrt(diag(xUncertainty)));

end

