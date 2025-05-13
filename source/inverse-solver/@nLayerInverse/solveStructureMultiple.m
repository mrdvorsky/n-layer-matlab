function [Parameters, Gamma, Uncertainty] = solveStructureMultiple(NLsolver, NL, f, gam, options)
%Perform simultaneous curve fitting on multiple nLayerInverse objects.
% This function takes quadruplets of nLayerInverse objects, nLayerForward
% objects, frequency vectors, and measurements, and tries to find the
% missing structure parameters of er, ur, thk to minimize the sum of the
% rms values of 'NL.calculate(f, er, ur, thk) - gam' for each set of
% quadruplets.
%
% This function is similar to the 'nLayerInverse.solveStructure' function
% except that each input set can have a different NLsolver object. See
% documentation of the 'solveStructure' function for more information.
% 
% Each NLsolver object must have the same number of parameters to be
% solved, but can otherwise be different (e.g., a different number of
% layers is ok). The values of er, ur, and thk passed into each NL object
% will be based on each corresponding NLsolver object, but the common
% parameters being solved for will be have the same value.
%
% As a filled-waveguide example, this function can be used to find a
% single permittivity for an MUT measured with two different thicknesses
% (or different waveguide bands, etc.). This is done by constructing two
% different nLayerInverse objects with the different MUT and waveguide
% section thicknesses, and setting the MUT erp and erpp to be solved for in
% both. Passing them into this function along with the measurements and NL
% objects will give the desired permittivity values.
%
% Example Usage (for multi-thickness, multi-band MUT measurement):
%   NLsolver1.setInitialValues(Er=[1, 4-0.1j, 1], Thk=[10, 2,  10]);
%   NLsolver2.setInitialValues(Er=[1, 4-0.1j, 1], Thk=[40, 10, 40]);
%   NLsolver1.setLayersToSolve(Er=[2]);
%   NLsolver2.setLayersToSolve(Er=[2]);
%   [Params, Gamma, Uncert] = nLayerInverse.solveStructureMultiple(...
%       NLsolver1, NL1, f1, gam1, ...
%       NLsolver2, NL2, f2, gam2);
%
% Example Usage (for multi-standoff open-ended measurements):
%   NLsolver1.setInitialValues(Er=[1, 2-0.01j], Thk=[0,  20]);
%   NLsolver2.setInitialValues(Er=[1, 2-0.01j], Thk=[10, 20]);
%   NLsolver1.setLayersToSolve(Er=[2], Thk=[2]);
%   NLsolver2.setLayersToSolve(Er=[2], Thk=[2]);
%   [Params, Gamma, Uncert] = nLayerInverse.solveStructureMultiple(...
%       NLsolver1, NL, f, gam1, ...
%       NLsolver2, NL, f, gam2);
%
%
% If the NLsolver objects have different min/max ranges and initial values
% for the common parameters, the initial values for the first NLsolver will
% be used, and the tightest min/max ranges for each parameter will be used.
% Additionally, the verbosity and local/global optimizer settings of the
% first NLsolver object will be used.
%
% Inputs:
%   NLsolver (Repeating) - A valid nLayerInverse object. Each must have
%       the same number of parameters to solve.
%   NL (Repeating) - A valid nLayerForward object.
%   f (Repeating) - Vector of frequencies to pass to NL.
%   gam (Repeating) - Measurements to fit. Size must match the output size
%       of 'NL.calculate(f, er, ur, thk)'.
%
% Outputs:
%   Parameters - Cell array of structs containing the structure parameters
%       (er, ur, thk), and simulated measurements (gam) for each input set.
%   Gamma - Cell array of simulated measurements (i.e., the value of
%       'NL.calculate(f, er, ur, thk)').
%   Uncertainty - Cell array of structs containing the calculated output
%       parmeter uncertainties for each input set.
%
% Named Arguments:
%   NoiseStdMin (0.01) - Minimum uncertainty value to assume for the
%       measurement data when calculating the Uncertainty struct. If the
%       RMS difference between the fit and the measurements is less than
%       NoiseStdMin, NoiseStdMin will be used instead.
%
% Author: Matt Dvorsky

arguments (Repeating)
    NLsolver(1, 1) nLayerInverse;
    NL(1, 1) nLayerForward;
    f(:, 1) {mustBeNonempty};
    gam {mustBeCorrectGamSize(f, gam)};
end
arguments
    options.NoiseStdMin(1, 1) {mustBeNonnegative} = 0.01;
end

%% Construct Linearized Ranges and Initial Guesses
[xInitial, xMin, xMax, xA, xb, xAeq, xbeq] = ...
    NLsolver{1}.constructInitialValuesAndRanges();
for ii = 2:numel(NLsolver)
    [~, xMinTmp, xMaxTmp] = ...
        NLsolver{ii}.constructInitialValuesAndRanges();
    xMin = max(xMin, xMinTmp);
    xMax = max(xMax, xMaxTmp);
end

%% Create Error Function
errorFunctionVector = @(x) nLayerInverse.calculateError(...
    NLsolver, x, NL, f, gam);
errorFunctionScalar = @(x) nLayerInverse.calculateError(...
    NLsolver, x, NL, f, gam, VectorOutput=false);

%% Set Verbosity for Optimizers
globalOptimizerOptions = NLsolver{1}.globalOptimizerOptions;
if NLsolver{1}.verbosity > 0
    globalOptimizerOptions = optimoptions(NLsolver{1}.globalOptimizerOptions, ...
        Display="iter");
end

localOptimizerOptions = NLsolver{1}.localOptimizerOptions;
if NLsolver{1}.verbosity > 0
    localOptimizerOptions = optimoptions(NLsolver{1}.localOptimizerOptions, ...
        Display="iter");
end

%% Run Global Optimizer
fit_res = 0;
if NLsolver{1}.useGlobalOptimizer && ~isempty(xInitial)
    if any(~isfinite(xMax))
        error("When the global optimizer option is enabled, " + ...
            "rangeMax_{er, ur, thk} must be set to finite values.");
    end
    switch class(globalOptimizerOptions)
        case "optim.options.GaOptions"
            [xInitial, fit_res] = ga(errorFunctionScalar, numel(xInitial), ...
                xA, xb, xAeq, xbeq, xMin, xMax, [], globalOptimizerOptions);
        case "optim.options.PatternsearchOptions"
            [xInitial, fit_res] = patternsearch(errorFunctionScalar, xInitial, ...
                xA, xb, xAeq, xbeq, xMin, xMax, [], globalOptimizerOptions);
        case "optim.options.Surrogateopt"
            [xInitial, fit_res] = surrogateopt(errorFunctionScalar, ...
                xMin, xMax, xA, xb, xAeq, xbeq, [], globalOptimizerOptions);
        case "optim.options.Particleswarm"
            if (numel(xA) + numel(xAeq)) > 0
                error("Thickness constraints not supported for the " + ...
                    "'Particleswarm' solver.");
            end
            [xInitial, fit_res] = particleswarm(errorFunctionScalar, numel(xInitial), ...
                xMin, xMax, globalOptimizerOptions);
        case "optim.options.SimulannealbndOptions"
            if (numel(xA) + numel(xAeq)) > 0
                error("Thickness constraints not supported for the " + ...
                    "'Particleswarm' solver.");
            end
            [xInitial, fit_res] = simulannealbnd(errorFunctionScalar, xInitial, ...
                xMin, xMax, globalOptimizerOptions);
        otherwise
            error("Global optimizer '%s' not supported.", ...
                class(globalOptimizerOptions));
    end
end

%% Run Local Optimizer
if NLsolver{1}.useLocalOptimizer && ~isempty(xInitial)
    switch class(localOptimizerOptions)
        case "optim.options.Lsqnonlin"
            [x, fit_res] = lsqnonlin_wrapper(errorFunctionVector, xInitial(:), ...
                xMin, xMax, xA, xb, xAeq, xbeq, localOptimizerOptions);
        case "optim.options.Fmincon"
            [x, fit_res] = fmincon(errorFunctionScalar, xInitial(:), ...
                xA, xb, xAeq, xbeq, xMin, xMax, [], localOptimizerOptions);
        otherwise
            error("Local optimizer '%s' not supported.", ...
                class(localOptimizerOptions));
    end
else
    x = xInitial(:);

    if isempty(xInitial)
        warning("No parameters were set to be optimized.");
    end
    
    if ~NLsolver{1}.useGlobalOptimizer
        warning("No local or global optimizer was enabled. " + ...
            "Results will be equal to initial conditions.");
    end
end

%% Check for Active Boundary Conditions
bounds_eps = 1e-6;
if any(abs(x - xMin) < bounds_eps, "all") || any(abs(x - xMax) < bounds_eps, "all")
    warning("One or more parameters (er, ur, thk) of the final " + ...
        "solved structure is bounded by a max or min, and thus the " + ...
        "solution is likely incorrect. Either relax the min and " + ...
        "max constraints or reduce the number of parameters " + ...
        "being solved.");
end

%% Create Parameters Output
Parameters = cell(numel(NLsolver), 1);
for ii = 1:numel(NLsolver)
    [er, ur, thk] = NLsolver{ii}.extractStructure(x, f);
    Parameters{ii}.er = er;
    Parameters{ii}.ur = ur;
    Parameters{ii}.thk = thk;
end

%% Create Gamma Output
if nargin >= 2
    Gamma = cell(numel(NLsolver), 1);
    for ii = 1:numel(NLsolver)
        Gamma{ii} = NL{ii}.calculate(f{ii}, Parameters{ii}.er, ...
            Parameters{ii}.ur, Parameters{ii}.thk);
    end
end

%% Create Uncertainty Output
if nargin >= 3
    inputParams = [NLsolver; NL; f];
    Uncertainty = nLayerInverse.computeParameterUncertaintyMultiple(...
        inputParams{:}, noiseStd=max(sqrt(fit_res), options.NoiseStdMin));
end

end




%% Helper Functions
function mustBeCorrectGamSize(f, gam)
    if iscell(f)    % Fix MATLAB bug.
        currentInd = find(cellfun(@(x) numel(x) > 0, f), 1, "last");
        f = f{currentInd};
    end
    if numel(f) ~= size(gam, 1)
        throwAsCaller(MException("nLayerInverse:mustBeCorrectGamSize", ...
            "First dimension of the measurements array must have size " + ...
            "equal to the number of frequencies (%d).", numel(f)));
    end
end

function [x, val] = lsqnonlin_wrapper(fun, x0, xMin, xMax, A, b, Aeq, beq, options)
    % Helper function to run lsqnonlin using Jacobian-based algorithm when
    % only equality constraints are specified. If inequality constraints
    % are specified or if there are no equality constraints, standard
    % lsqnonlin will be used. Otherwise, a change of variables will be used
    % to reduce the number of unknowns. Unfortunately, MATLAB currently has
    % no algorithm that does this when using a Jacobian-based solver.
    
    % Run normally if there are any inequality or no equality constraints.
    if (numel(A) > 0) || (numel(Aeq) == 0)
        [x, val] = lsqnonlin(fun, x0, xMin, xMax, ...
            A, b, Aeq, beq, [], options);
        return;
    end
    
    % If only equality constraints are specified, use substitution method.
    [R, xSubInd(:, 1)] = rref([Aeq, beq]);
    yInd(:, 1) = setdiff((1:numel(x0)).', xSubInd);

    % Write x in the form of C*y + c0.
    C = zeros(numel(x0), numel(yInd));
    C(yInd, :) = eye(numel(yInd));
    C(xSubInd, :) = -R(:, yInd);

    c0 = zeros(numel(x0), 1);
    c0(xSubInd, 1) = R(:, end);

    % Identify y0, yMin, and yMax.
    yMin = xMin(yInd);
    yMax = xMax(yInd);
    y0 = C \ (x0 - c0);

    [y, val] = lsqnonlin(@(y) fun(max(min(C*y + c0, xMax), xMin)), ...
        y0, yMin, yMax, [], [], [], [], [], options);

    x = C*y + c0;
end

