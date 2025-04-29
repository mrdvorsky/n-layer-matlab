classdef nLayerInverse < matlab.mixin.Copyable
    %NLAYERINVERSE Class to perform optimization on "nLayerForward" objects.
    % This class allows a user to define an nLayer structure on which to
    % perform curve fitting.
    %
    % Example Usage: (check examples folder for more use cases)
    %   NL = ...;           % Forward solver. Can be any "nLayerForward".
    %
    %   NLsolver = nLayerInverse(2, verbosity=1);   % Create 2-layer structure.
    %   NLsolver.setInitialValues(Er=erStart, Thk=thkStart);
    %   NLsolver.setLayersToSolve(Er=[1], Thk=[1, 2]);
    %
    %   NLsolver.printStructureParameters(showLimits=true);
    %   [Params, ~, Uncert] = ...
    %       NLsolver.solveStructure(NL, f, gam);    % Minimize "NL.calculate(f, ...) - gam"
    %   NLsolver.printStructureParameters(Params, Uncert);
    %
    % nLayerInverse Properties:
    %   verbosity (0) - Verbosity level. Set to 1 to enable the optimizer
    %       console output.
    %
    % Author: Matt Dvorsky

    properties (SetAccess=immutable)
        layerCount(1, 1) {mustBeInteger, mustBePositive} = 1;   % Number of layers in the structure.
    end
    properties (Access=public)
        initialValue_er(1, :) {} = 1;               % Structure er values.
        initialValue_ur(1, :) {} = 1;               % Structure ur values.
        initialValue_thk(1, :) {mustBeNonnegative} = 0;     % Structure thk values.

        layersToSolve_erp(1, :)  {mustBeInteger, mustBePositive} = [];  % Layer indices to solve for erp.
        layersToSolve_erpp(1, :) {mustBeInteger, mustBePositive} = [];  % Layer indices to solve for erpp.
        layersToSolve_urp(1, :)  {mustBeInteger, mustBePositive} = [];  % Layer indices to solve for urp.
        layersToSolve_urpp(1, :) {mustBeInteger, mustBePositive} = [];  % Layer indices to solve for urpp.
        layersToSolve_thk(1, :)  {mustBeInteger, mustBePositive} = [];  % Layer indices to solve for thk.

        rangeMin_erp(1, :)  {mustBeReal} = 1;       % Minimum erp per layer.
        rangeMin_erpp(1, :) {mustBeReal} = 0;       % Minimum erpp per layer.
        rangeMin_urp(1, :)  {mustBeReal} = 1;       % Minimum urp per layer.
        rangeMin_urpp(1, :) {mustBeReal} = 0;       % Minimum urpp per layer.
        rangeMin_thk(1, :)  {mustBeReal} = 0;       % Minimum thk per layer.

        rangeMax_erp(1, :)  {mustBeReal} = inf;     % Maximum erp per layer.
        rangeMax_erpp(1, :) {mustBeReal} = inf;     % Maximum erpp per layer.
        rangeMax_urp(1, :)  {mustBeReal} = inf;     % Maximum urp per layer.
        rangeMax_urpp(1, :) {mustBeReal} = inf;     % Maximum urpp per layer.
        rangeMax_thk(1, :)  {mustBeReal} = inf;     % Maximum thk per layer.

        verbosity(1, 1) {mustBeInteger, mustBeNonnegative} = 0; % Verbosity level. Set to 0 for no output.

        useGlobalOptimizer(1, 1) logical = false;   % Whether or not to use global optimizer.
        useLocalOptimizer(1, 1) logical = true;     % Whether or not to use local optimizer.

        localOptimizerOptions(1, 1) ...                         % Options object for local optimizer.
            {mustBeA(localOptimizerOptions, "optim.options.SolverOptions")} ...
            = optimoptions("lsqnonlin", Display="none");
        globalOptimizerOptions(1, 1) ...                        % Options object for global optimizer.
            {mustBeA(globalOptimizerOptions, "optim.options.SolverOptions")} ...
            = optimoptions("surrogateopt", Display="none", PlotFcn="");
    end
    properties (GetAccess=public, SetAccess=protected)
        constraints_thk_Aeq(:, :) = [];
        constraints_thk_A(:, :) = [];
        constraints_thk_b(:, 1) = [];
    end

    %% Class Constructor
    methods
        function O = nLayerInverse(numLayers, classProperties)
            %NLAYERINVERSE Construct an instance of this class.
            %
            % Inputs:
            %   numLayers - Number of structure layers.

            arguments
                numLayers(1, 1) {mustBeInteger, mustBePositive};
            end
            arguments
                classProperties.?nLayerInverse;
            end

            % Set number of layers.
            O.layerCount = numLayers;
            O.initialValue_er = repmat(O.initialValue_er, numLayers, 1);
            O.initialValue_ur = repmat(O.initialValue_ur, numLayers, 1);
            O.initialValue_thk = repmat(O.initialValue_thk, numLayers, 1);
            O.setRanges();
            
            O.constraints_thk_Aeq = zeros(0, numLayers);
            O.constraints_thk_A = zeros(0, numLayers);
            O.constraints_thk_b = zeros(0, 1);

            % Set class parameter values.
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                O.(propPairs{ii}) = propPairs{ii + 1};
            end
        end
    end

    %% Class Functions
    methods (Access=public)
        setLayersToSolve(O, options);
        setInitialValues(O, options);
        setRanges(O, options);
        addThicknessConstraint(O, layerIndices, constraints);

        [varargout] = solveStructure(O, NL, f, gam, options);
        [Uncert] = computeParameterUncertainty(O, NL, f, options);

        [varargout] = printStructureParameters(O, Parameters, Uncertainty, formatOptions, options);
    end
    methods (Access=private)
        [xGuess, xMin, xMax, xA, xb, xAeq, xbeq] = ...
            constructInitialValuesAndRanges(O);
        [er, ur, thk] = extractStructure(O, x, f);
        validate(O);
    end
    methods (Static, Access=public)
        [Params, Gamma, Uncert] = solveStructureMultiple(NLsolver, NL, f, gam, options);
        [Uncert] = computeParameterUncertaintyMultiple(NLsolver, NL, f, options);
    end
    methods (Static, Access=private)
        [gamError, gamErrorComplex] = calculateError(...
            NLsolver, x, NL, f, gamActual, options);
    end

    %% Class Setters
    methods
        % Setters for initial structure values.
        function set.initialValue_er(O, er)
            checkInitialValues(O, er);
            O.initialValue_er = er;
        end
        function set.initialValue_ur(O, ur)
            checkInitialValues(O, ur);
            O.initialValue_ur = ur;
        end
        function set.initialValue_thk(O, thk)
            checkInitialValues(O, thk);
            O.initialValue_thk = thk;
        end

        % Setters for layers to solve.
        function set.layersToSolve_erp(O, erp)
            checkLayersToSolve(O, erp);
            O.layersToSolve_erp = sort(erp);
        end
        function set.layersToSolve_erpp(O, erpp)
            checkLayersToSolve(O, erpp);
            O.layersToSolve_erpp = sort(erpp);
        end
        function set.layersToSolve_urp(O, urp)
            checkLayersToSolve(O, urp);
            O.layersToSolve_urp = sort(urp);
        end
        function set.layersToSolve_urpp(O, urpp)
            checkLayersToSolve(O, urpp);
            O.layersToSolve_urpp = sort(urpp);
        end
        function set.layersToSolve_thk(O, thk)
            checkLayersToSolve(O, thk);
            O.layersToSolve_thk = sort(thk);
        end

        % Setters for minimum values of range.
        function set.rangeMin_erp(O, erp)
            erp = checkRanges(O, erp);
            O.rangeMin_erp = erp;
        end
        function set.rangeMin_erpp(O, erpp)
            erpp = checkRanges(O, erpp);
            O.rangeMin_erpp = erpp;
        end
        function set.rangeMin_urp(O, urp)
            urp = checkRanges(O, urp);
            O.rangeMin_urp = urp;
        end
        function set.rangeMin_urpp(O, urpp)
            urpp = checkRanges(O, urpp);
            O.rangeMin_urpp = urpp;
        end
        function set.rangeMin_thk(O, thk)
            thk = checkRanges(O, thk);
            O.rangeMin_thk = thk;
        end

        % Setters for maximum values of range.
        function set.rangeMax_erp(O, erp)
            erp = checkRanges(O, erp);
            O.rangeMax_erp = erp;
        end
        function set.rangeMax_erpp(O, erpp)
            erpp = checkRanges(O, erpp);
            O.rangeMax_erpp = erpp;
        end
        function set.rangeMax_urp(O, urp)
            urp = checkRanges(O, urp);
            O.rangeMax_urp = urp;
        end
        function set.rangeMax_urpp(O, urpp)
            urpp = checkRanges(O, urpp);
            O.rangeMax_urpp = urpp;
        end
        function set.rangeMax_thk(O, thk)
            thk = checkRanges(O, thk);
            O.rangeMax_thk = thk;
        end
    end
end




%% Helper Functions
function checkInitialValues(O, erUrThk)
    if O.layerCount ~= numel(erUrThk)
        error("The parameters 'initialValue_{er, ur, thk}' " + ...
            "must be vectors with layerCount (%d) elements.", ...
            O.layerCount);
    end
end

function checkLayersToSolve(O, inds)
    if any(inds > O.layerCount, "all") || any(inds <= 0, "all") ...
            || numel(inds) ~= numel(unique(inds)) ...
            || any(inds ~= floor(inds), "all")
        error("The parameters 'layersToSolve_{er, ur, thk}' must " + ...
            "consist of unique positive integers no greater than " + ...
            "the layer count (%d).", O.layerCount);
    end
end

function [rangeVal] = checkRanges(O, rangeVal)
    if (O.layerCount ~= numel(rangeVal)) && (1 ~= numel(rangeVal))
        error("The parameters 'range{Min, Max}_{er, ur, thk}' must " + ...
            "be scalars or vectors with layerCount (%d) elements.", ...
            O.layerCount);
    end
    if numel(rangeVal) == 1
        rangeVal = repmat(rangeVal, O.layerCount, 1);
    end
end

