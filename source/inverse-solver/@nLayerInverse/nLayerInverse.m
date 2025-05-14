classdef nLayerInverse < matlab.mixin.Copyable
    %Class to perform optimization on "nLayerForward" objects.
    % This class allows a user to define an nLayer structure on which to
    % perform curve fitting.
    %
    % ===== Basic Usage =====
    %   NL = ...;       % Forward solver. Can be any "nLayerForward".
    %
    %   NLsolver = nLayerInverse(2, verbosity=1);   % 2-layer structure.
    %   NLsolver.setInitialValues(Er=erStart, Thk=thkStart);
    %   NLsolver.setLayersToSolve(Er=[1], Thk=[1, 2]);
    %
    %   NLsolver.printStructureParameters(showLimits=true);
    %   
    %   % Minimize "NL.calculate(f, ...) - gam"
    %   [Params, ~, Uncert] = NLsolver.solveStructure(NL, f, gam);
    %
    %   NLsolver.printStructureParameters(Params, Uncert);
    %
    %
    % Author: Matt Dvorsky

    properties (SetAccess=immutable)
        layerCount(1, 1) {mustBeInteger, mustBePositive} = 1;   % Number of layers in the structure.
    end
    properties (Access=public)
        initialValue_er(1, :) {} = 1;               % Structure er values.
        initialValue_ur(1, :) {} = 1;               % Structure ur values.
        initialValue_thk(1, :) {} = 0;     % Structure thk values.

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
        function self = nLayerInverse(numLayers, classProperties)
            %Construct an instance of this class.
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
            self.layerCount = numLayers;
            self.initialValue_er = repmat(self.initialValue_er, numLayers, 1);
            self.initialValue_ur = repmat(self.initialValue_ur, numLayers, 1);
            self.initialValue_thk = repmat(self.initialValue_thk, numLayers, 1);
            self.setRanges();
            
            self.constraints_thk_Aeq = zeros(0, numLayers);
            self.constraints_thk_A = zeros(0, numLayers);
            self.constraints_thk_b = zeros(0, 1);

            % Set class parameter values.
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                self.(propPairs{ii}) = propPairs{ii + 1};
            end
        end
    end

    %% Class Functions
    methods (Access=public)
        setLayersToSolve(self, options);
        setInitialValues(self, options);
        setRanges(self, options);
        addThicknessConstraint(self, layerIndices, constraints);

        [varargout] = solveStructure(self, NL, f, gam, options);
        [Uncert] = computeParameterUncertainty(self, NL, f, options);

        [varargout] = printStructureParameters(self, Parameters, ...
            Uncertainty, formatOptions, options);
    end
    methods (Access=private)
        [xGuess, xMin, xMax, xA, xb, xAeq, xbeq] = ...
            constructInitialValuesAndRanges(self);
        [er, ur, thk] = extractStructure(self, x, f);
        validate(self);
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
        function set.initialValue_er(self, er)
            er = checkInitialValues(self, er);
            self.initialValue_er = er;
        end
        function set.initialValue_ur(self, ur)
            ur = checkInitialValues(self, ur);
            self.initialValue_ur = ur;
        end
        function set.initialValue_thk(self, thk)
            thk = checkInitialValues(self, thk);
            self.initialValue_thk = thk;
            mustBeNonnegative(thk);
        end

        % Setters for layers to solve.
        function set.layersToSolve_erp(self, erp)
            checkLayersToSolve(self, erp);
            self.layersToSolve_erp = sort(erp);
        end
        function set.layersToSolve_erpp(self, erpp)
            checkLayersToSolve(self, erpp);
            self.layersToSolve_erpp = sort(erpp);
        end
        function set.layersToSolve_urp(self, urp)
            checkLayersToSolve(self, urp);
            self.layersToSolve_urp = sort(urp);
        end
        function set.layersToSolve_urpp(self, urpp)
            checkLayersToSolve(self, urpp);
            self.layersToSolve_urpp = sort(urpp);
        end
        function set.layersToSolve_thk(self, thk)
            checkLayersToSolve(self, thk);
            self.layersToSolve_thk = sort(thk);
        end

        % Setters for minimum values of range.
        function set.rangeMin_erp(self, erp)
            erp = checkRanges(self, erp);
            self.rangeMin_erp = erp;
        end
        function set.rangeMin_erpp(self, erpp)
            erpp = checkRanges(self, erpp);
            self.rangeMin_erpp = erpp;
        end
        function set.rangeMin_urp(self, urp)
            urp = checkRanges(self, urp);
            self.rangeMin_urp = urp;
        end
        function set.rangeMin_urpp(self, urpp)
            urpp = checkRanges(self, urpp);
            self.rangeMin_urpp = urpp;
        end
        function set.rangeMin_thk(self, thk)
            thk = checkRanges(self, thk);
            self.rangeMin_thk = thk;
        end

        % Setters for maximum values of range.
        function set.rangeMax_erp(self, erp)
            erp = checkRanges(self, erp);
            self.rangeMax_erp = erp;
        end
        function set.rangeMax_erpp(self, erpp)
            erpp = checkRanges(self, erpp);
            self.rangeMax_erpp = erpp;
        end
        function set.rangeMax_urp(self, urp)
            urp = checkRanges(self, urp);
            self.rangeMax_urp = urp;
        end
        function set.rangeMax_urpp(self, urpp)
            urpp = checkRanges(self, urpp);
            self.rangeMax_urpp = urpp;
        end
        function set.rangeMax_thk(self, thk)
            thk = checkRanges(self, thk);
            self.rangeMax_thk = thk;
        end
    end
end




%% Helper Functions
function [erUrThk] = checkInitialValues(self, erUrThk)
    if iscell(erUrThk)
        if ~all(cellfun(@isscalar, erUrThk))
            error("nLayerInverse:cellArrayWithNonScalars", ...
                "If the parameter 'initialValue_{er, ur, thk}' " + ...
                "is a cell array, it must be a vector of scalars.");
        end
        erUrThk = cell2mat(erUrThk);
    end
    if self.layerCount ~= numel(erUrThk)
        error("nLayerInverse:layerCountMismatch", ...
            "The parameters 'initialValue_{er, ur, thk}' " + ...
            "must be vectors with layerCount (%d) elements.", ...
            self.layerCount);
    end
end

function checkLayersToSolve(self, inds)
    if any(inds > self.layerCount, "all") || any(inds <= 0, "all") ...
            || numel(inds) ~= numel(unique(inds)) ...
            || any(inds ~= floor(inds), "all")
        error("nLayerInverse:layersToSolveImproperIndices", ...
            "The parameters 'layersToSolve_{er, ur, thk}' must " + ...
            "consist of unique positive integers no greater than " + ...
            "the layer count (%d).", self.layerCount);
    end
end

function [rangeVal] = checkRanges(self, rangeVal)
    if (self.layerCount ~= numel(rangeVal)) && (1 ~= numel(rangeVal))
        error("nLayerInverse:layerCountMismatch", ...
            "The parameters 'range{Min, Max}_{er, ur, thk}' must " + ...
            "be scalars or vectors with layerCount (%d) elements.", ...
            self.layerCount);
    end
    if isscalar(rangeVal)
        rangeVal = repmat(rangeVal, self.layerCount, 1);
    end
end

