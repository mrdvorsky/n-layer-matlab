function [er, ur, thk, varargout] = solveStructure(O, NL, f, gam)
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

%% Set Verbosity for Optimizers
globalOptimizerOptions = O.globalOptimizerOptions;
if O.verbosity > 0
    globalOptimizerOptions = optimoptions(O.globalOptimizerOptions, ...
        Display="iter");
end

localOptimizerOptions = O.localOptimizerOptions;
if O.verbosity > 0
    localOptimizerOptions = optimoptions(O.localOptimizerOptions, ...
        Display="iter");
end

%% Run Global Optimizer
if O.useGlobalOptimizer
    switch class(globalOptimizerOptions)
        case "optim.options.GaOptions"
            xInitial = ga(errorFunctionScalar, numel(xInitial), ...
                [], [], [], [], xMin, xMax, [], globalOptimizerOptions);
        case "optim.options.Particleswarm"
            xInitial = particleswarm(errorFunctionScalar, numel(xInitial), ...
                xMin, xMax, globalOptimizerOptions);
        case "optim.options.PatternsearchOptions"
            xInitial = patternsearch(errorFunctionScalar, xInitial, ...
                [], [], [], [], xMin, xMax, [], globalOptimizerOptions);
        case "optim.options.SimulannealbndOptions"
            xInitial = simulannealbnd(errorFunctionScalar, xInitial, ...
                xMin, xMax, globalOptimizerOptions);
        case "optim.options.Surrogateopt"
            xInitial = surrogateopt(errorFunctionScalar, ...
                xMin, xMax, [], [], [], [], [], globalOptimizerOptions);
        otherwise
            error("Global optimizer '%s' not supported.", ...
                class(globalOptimizerOptions));
    end
end

%% Run Local Optimizer
if O.useLocalOptimizer
    switch class(localOptimizerOptions)
        case "optim.options.Lsqnonlin"
            x = lsqnonlin(errorFunctionVector, xInitial, xMin, xMax, ...
                localOptimizerOptions);
        case "optim.options.Fmincon"
            x = fmincon(errorFunctionScalar, xInitial, ...
                [], [], [], [], xMin, xMax, [], localOptimizerOptions);
        otherwise
            error("Local optimizer '%s' not supported.", ...
                class(localOptimizerOptions));
    end
else
    x = xInitial;
    
    if ~O.useGlobalOptimizer
        error("At least one optimizer must be enabled.");
    end
end

%% Create Output
[er, ur, thk] = O.extractStructure(x, f);

%% Assign Gamma
if nargout > 3
    varargout{1} = NL.calculate(f, er, ur, thk);
end

end
