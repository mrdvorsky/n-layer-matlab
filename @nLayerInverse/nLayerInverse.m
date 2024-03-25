classdef nLayerInverse < matlab.mixin.Copyable & matlab.mixin.SetGetExactNames
    %NLAYERINVERSE Class to perform optimization on nLayerForward objects.
    % This class allows a user to define an nLayer structure on which to
    % perform curve fitting.
    %
    % Example Usage:
    %   fas
    %
    % nLayerInverse Properties:
    %   verbosity (0) - Verbosity level. Set to 1 to enable the optimizer
    %       console output.
    %
    % Author: Matt Dvorsky
    
    properties (Access=public)
        layersToSolve_erp(:, 1)  {mustBeInteger, mustBePositive} = [];  % Layer indices to solve for erp.
        layersToSolve_erpp(:, 1) {mustBeInteger, mustBePositive} = [];  % Layer indices to solve for erpp.
        layersToSolve_urp(:, 1)  {mustBeInteger, mustBePositive} = [];  % Layer indices to solve for urp.
        layersToSolve_urpp(:, 1) {mustBeInteger, mustBePositive} = [];  % Layer indices to solve for urpp.
        layersToSolve_thk(:, 1)  {mustBeInteger, mustBePositive} = [];  % Layer indices to solve for thk.
        
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

        initialValue_er(1, :) {} = 1;               % Structure er values.
        initialValue_ur(1, :) {} = 1;               % Structure ur values.
        initialValue_thk(1, :) {mustBeReal} = 0;    % Structure thk values.

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
    properties (GetAccess=public, SetAccess=private)
        layerCount(1, 1) {mustBeInteger, mustBePositive} = 1;   % Number of layers in the structure.
    end
    
    %% Class Functions
    methods (Access=public)
        setLayerCount(O, layerCount);
        setLayersToSolve(O, options);
        setInitialValues(O, options);
        setRanges(O, options);

        [varargout] = solveStructure(O, NL, f, gam, options);

        [Uncert] = computeParameterUncertainty(O, NL, f, options);
        [varargout] = printStructureParameters(O, Parameters, Uncertainty, formatOptions, options);
        validate(O);
    end
    methods (Access=private)
        [xGuess, xMin, xMax] = constructInitialValuesAndRanges(O);
        [er, ur, thk] = extractStructure(O, x, f);
    end
    methods (Static, Access=public)
        [Params, Gamma, Uncert] = solveStructureMultiple(NLsolver, NL, f, gam, options);
        [Uncert] = computeParameterUncertaintyMultiple(NLsolver, NL, f, options);
    end
    methods (Static, Access=private)
        [gamError, gamErrorComplex] = calculateError(...
            NLsolver, x, NL, f, gamActual, options);
    end
    
    %% Class Constructor
    methods
        function O = nLayerInverse(numLayers, classProperties)
            %NLAYERINVERSE Construct an instance of this class.
            % Example Usage:
            %   See example usage in main class documentation. Note that
            %   all public class properties can be specified as a named
            %   argument to the constructor (e.g., as "verbosity=1").
            %
            % Inputs:
            %   numLayers - Number of structure layers.
            
            arguments
                numLayers(1, 1) {mustBeInteger, mustBePositive};
            end
            arguments (Repeating)
                classProperties;
            end

            % Set Class Parameter Values
            if mod(numel(classProperties), 2) ~= 0
                error("Parameter and value arguments must come in pairs.");
            end
            for ii = 1:2:numel(classProperties)
                set(O, classProperties{ii}, classProperties{ii + 1});
            end
            
            % Set Number of Layers
            O.setLayerCount(numLayers);
        end
    end
    
end

