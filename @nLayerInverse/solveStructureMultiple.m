function [Parameters, Gamma, Uncertainty] = solveStructureMultiple(NLsolver, NL, f, gam)
%SOLVESTRUCTUREMULTI Perform simultaneous curve fitting on multiple nLayerInverse objects.
% This function takes quadruplets of nLayerInverse objects, nLayerForward
% objects, frequency vectors, and measurements, and tries to find the
% missing structure parameters of er, ur, thk to minimize the sum of the
% rms values of 'NL.calculate(er, ur, thk) - gam' for each quadruplets.
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
%   NLsolver1.setLayersToSolve(Erp=[2], Erpp=[2]);
%   NLsolver2.setLayersToSolve(Erp=[2], Erpp=[2]);
%   [Params, Uncert] = nLayerInverse.solveStructureMultiple(...
%       NLsolver1, NL1, f1, gam1, ...
%       NLsolver2, NL2, f2, gam2);
%
% Example Usage (for multi-standoff open-ended measurements):
%   NLsolver1.setInitialValues(Er=[1, 2-0.01j], Thk=[0,  20]);
%   NLsolver2.setInitialValues(Er=[1, 2-0.01j], Thk=[10, 20]);
%   NLsolver1.setLayersToSolve(Erp=[2], Erpp=[2], Thk=[2]);
%   NLsolver2.setLayersToSolve(Erp=[2], Erpp=[2], Thk=[2]);
%   [Params, Uncert] = nLayerInverse.solveStructureMultiple(...
%       NLsolver1, NL, f, gam1, ...
%       NLsolver2, NL, f, gam2);
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
% Author: Matt Dvorsky

arguments (Repeating)
    NLsolver(1, 1) {mustBeA(NLsolver, "nLayerInverse")};
    NL(1, 1) {mustBeA(NL, "nLayerForward")};
    f(:, 1) {mustBeNonempty};
    gam(:, :) {mustBeCorrectGamSize(f, gam)};
end

%% Validate nLayerInverse Objects
for ii = 1:numel(NLsolver)
    NLsolver{ii}.validate();
end

%% Construct Linearized Ranges and Initial Guesses
[xInitial, xMin, xMax] = NLsolver{1}.constructInitialValuesAndRanges();
for ii = 2:numel(NLsolver)
    [~, xMinTmp, xMaxTmp] = NLsolver{ii}.constructInitialValuesAndRanges();
    xMin = max(xMin, xMinTmp);
    xMax = max(xMax, xMaxTmp);
end

%% Create Error Function
errorFunctionVector = @(x) nLayerInverse.calculateError(NLsolver, x, NL, f, gam);
errorFunctionScalar = @(x) nLayerInverse.calculateError(NLsolver, x, NL, f, gam, ...
    VectorOutput=false);

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
if NLsolver{1}.useGlobalOptimizer
    if any(~isfinite(xMax))
        error("When the global optimizer option is enabled, " + ...
            "rangeMax_{er, ur, thk} must be set to finite values.");
    end
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
if NLsolver{1}.useLocalOptimizer
    switch class(localOptimizerOptions)
        case "optim.options.Lsqnonlin"
            x = lsqnonlin(errorFunctionVector, xInitial(:), xMin, xMax, ...
                localOptimizerOptions);
        case "optim.options.Fmincon"
            x = fmincon(errorFunctionScalar, xInitial(:), ...
                [], [], [], [], xMin, xMax, [], localOptimizerOptions);
        otherwise
            error("Local optimizer '%s' not supported.", ...
                class(localOptimizerOptions));
    end
else
    x = xInitial(:);
    
    if ~NLsolver{1}.useGlobalOptimizer
        error("At least one optimizer must be enabled.");
    end
end

%% Check for Active Boundary Conditions
bounds_eps = 1e-8;
if any(abs(x - xMin) < bounds_eps) || any(abs(x - xMax) < bounds_eps)
    warning("One or more parameters (er, ur, thk) of the final " + ...
        "solved structure is bounded by a max or min, and thus the " + ...
        "solution is likely incorrect. Either relax the min and " + ...
        "max constraints or reduce the number of parameters " + ...
        "being solved.")
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
    Uncertainty = struct();
end

end


function mustBeCorrectGamSize(f, gam)
    if iscell(f)    % Fix MATLAB bug.
        currentInd = find(cellfun(@(x) numel(x) > 0, f), 1, "last");
        f = f{currentInd};
    end
    if numel(f) ~= size(gam)
        throwAsCaller(MException("nLayerInverse:mustBeCorrectGamSize", ...
            "First dimension of the measurements array must have size " + ...
            "equal to the number of frequencies (%d).", numel(f)));
    end
end

