classdef nLayerInverse < matlab.mixin.Copyable
    %NLAYERINVERSE Class to perform optimization on nLayerForward objects.
    % This class serves as an interface definition for all nLayer forward
    % calculator objects. These objects take in a multilayer structure
    % definition and output the computed reflection coefficients (or
    % transmission coefficients). It also contains several useful utility
    % functions.
    %
    % nLayerForward Properties:
    %   speedOfLight (299.792458) - Speed of light (mm GHz). Must match
    %       units of distance and frequency used.
    %   verbosity (0) - A value of 0 should suppress console output.
    %   checkStructureValues (true) - Flag used in the "verifyStructure"
    %       function. If true, this function will throw errors if
    %       non-physical values of er, ur, or thk are passed in.
    %
    % Author: Matt Dvorsky
    
    properties (GetAccess = public, SetAccess = public)
        localOptimizerOptions;
        globalOptimizerOptions;
        
        useGlobalOptimizer = false;
        useLocalOptimizer = true;
        
        default_erRange = [1; 10];
        default_erpRange = [0.001; 10];
        default_thkRange = [0.001; 1];
    end
    
    properties (GetAccess = public, SetAccess = private)
        verbosity;
        
        layerCount;
        
        erLayersToSolve;
        erpLayersToSolve;
        thkLayersToSolve;
        
        erRange;
        erpRange;
        thkRange;
        
        erInitialValue;
        erpInitialValue;
        thkInitialValue;
    end
    
    %% Function Declarations (implemented in separate files)
    methods (Access = public)
        [er, ur, thk, varargout] = solveStructure(O, NL, f, gam);
        setLayerCount(O, layerCount);
        setLayersToSolve(O, options);
        setInitialValues(O, options);
        [varargout] = printStructureParameters(O, options);
    end
    
    methods (Access = private)
        [er, ur, thk] = extractStructure(O, x, f);
        [xGuess, xMin, xMax] = constructInitialValuesAndRanges(O);
    end
    
    %% Class constructor
    methods
        function O = nLayerInverse(numLayers, options)
            %NLAYERRECTANGULAR Construct an instance of this class.
            % Example Usage:
            %   NLsolver = nLayerInverse();
            %
            % Inputs:
            %   maxM - .
            % Named Arguments:
            %   LocalOptimizerOptions - .
            
            arguments
                numLayers(1, 1) {mustBeInteger, mustBePositive};
                options.Verbosity = 0;
                options.LocalOptimizerOptions;
                options.GlobalOptimizerOptions;
            end
            
            %% Set Class Parameter Values
            O.verbosity = options.Verbosity;
            
            %% Set Default Optimizer Options
            if isfield(options, "LocalOptimizerOptions")
                O.localOptimizerOptions = options.LocalOptimizerOptions;
            else
                O.localOptimizerOptions = optimoptions("lsqnonlin", ...
                    Display="none");
            end
            
            if isfield(options, "GlobalOptimizerOptions")
                O.globalOptimizerOptions = options.GlobalOptimizerOptions;
            else
                O.globalOptimizerOptions = optimoptions("surrogateopt", ...
                    Display="none", PlotFcn="");
            end
            
            %% Set Number of Layers
            O.setLayerCount(numLayers);
            
        end
    end
    
end

