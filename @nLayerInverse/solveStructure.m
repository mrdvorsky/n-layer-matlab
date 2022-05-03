function [er, ur, thk] = solveStructure(O, NL, f, gam)
%SOLVESTRUCTURE Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
    NL;
    f(:, 1);
    gam;
end

%% Construct Linearized Ranges and Initial Guesses
[xInitial, xMin, xMax] = O.constructInitialValuesAndRanges();

%% Create Error Function
errorFunctionVector = @(x) O.calculateError(x, NL, f, gam);
errorFunctionScalar = @(x) O.calculateError(x, NL, f, gam, ...
    VectorOutput=false);

%% Run Global Optimizer
if O.useGlobalOptimizer
    switch class(O.globalOptimizerOptions)
        case "optim.options.GaOptions"
            xInitial = ga(errorFunctionScalar, numel(xInitial), ...
                [], [], [], [], xMin, xMax, [], O.globalOptimizerOptions);
        case "optim.options.Particleswarm"
            xInitial = particleswarm(errorFunctionScalar, numel(xInitial), ...
                xMin, xMax, O.globalOptimizerOptions);
        case "optim.options.PatternsearchOptions"
            xInitial = patternsearch(errorFunctionScalar, xInitial, ...
                [], [], [], [], xMin, xMax, [], O.globalOptimizerOptions);
        case "optim.options.SimulannealbndOptions"
            xInitial = simulannealbnd(errorFunctionScalar, xInitial, ...
                xMin, xMax, O.globalOptimizerOptions);
        case "optim.options.Surrogateopt"
            xInitial = ga(errorFunctionScalar, numel(xInitial), ...
                [], [], [], [], xMin, xMax, [], O.globalOptimizerOptions);
        otherwise
            error("Global optimizer '%s' not supported.", ...
                class(O.globalOptimizerOptions));
    end
end

%% Run Local Optimizer
if O.useLocalOptimizer
    switch class(O.localOptimizerOptions)
        case "optim.options.Lsqnonlin"
            x = lsqnonlin(errorFunctionVector, xInitial, xMin, xMax, ...
                O.localOptimizerOptions);
        case "optim.options.Fmincon"
            x = fmincon(errorFunctionScalar, xInitial, ...
                [], [], [], [], xMin, xMax, [], O.localOptimizerOptions);
        otherwise
            error("Local optimizer '%s' not supported.", ...
                class(O.localOptimizerOptions));
    end
else
    x = xInitial;
    
    if ~O.useGlobalOptimizer
        error("At least one optimizer must be enabled.");
    end
end

%% Create Output
[er, ur, thk] = O.extractStructure(x, f);

end
