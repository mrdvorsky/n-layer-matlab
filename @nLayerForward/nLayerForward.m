classdef (Abstract) nLayerForward < matlab.mixin.Copyable
    %NLAYERFORWARD Interface class for nLayer forward calculators.
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
    % To use this class, subclass it. Any subclasses must, at a minimum,
    % implement the "calculateGamma" function and the "getOutputLabels"
    % function, whose functionalities are described below.
    % 
    % Abstract functions (should be implemented by subclass):
    %   calculateGamma - Should accept a vector of frequencies and the
    %       multilayer structure definition at each frequency, and should
    %       output the calculated reflection (and/or transmission)
    %       coefficient(s) at each frequency. The first column of the
    %       output should have the same length as the number of
    %       frequencies. Additional dimensions can be anything.
    %   getOutputLabels - Should return a vector of strings labeling each
    %       output channel. If "calculateGamma" returns gam, then the nth
    %       element of this output vector should describe gam(:, n).
    %
    % The implementation of "calculateGamma" should utilize the parameter
    % "speedOfLight" for defining units (default mm GHz) and observe the
    % "verbosity" parameter value.
    %
    % Utility functions (implemented by this class):
    %   fejer2 - Generates weights and nodes to perform Fejer Type II
    %       quadrature integration.
    %   gaussKronrod - Generates weights and nodes to perform Gaussian
    %       quadrature integration.
    %   integralVectorized - Routine to quickly perform adaptive integration
    %       of vectorized functions.
    %   verifyStructure - Checks the validity of the multilayer structure
    %       and frequency definitions. This function is called
    %       automatically on inputs passed into calculate(...).
    %   printStructure - Prints or returns a string to visualize the
    %       multilayer structure.
    %   changeStructureConductivity - Changes a multilayer structure to be
    %       conductor-backed with a finite conductivity.
    %
    % Author: Matt Dvorsky
    
    properties (GetAccess = public, SetAccess = public)
        speedOfLight = 299.792458;      % Speed of light (mm/ns).
        verbosity = 0;                  % Verbosity level. Zero for no console output.
        checkStructureValues = true;    % Whether to check ranges of er, ur, and thk.
    end
    
    %% Virtual Protected member function definitions
    methods (Abstract, Access = protected)
        [gam] = calculateGamma(O, f, er, ur, thk);
    end
    
    %% Virtual Public member function definitions
    methods (Abstract, Access = public)
        [outputLabels] = getOutputLabels(O);
    end
    
    %% Public member function definitions (implemented in separate files)
    methods (Access = public)
        [gam] = calculate(O, f, er, ur, thk);
        [er, ur, thk] = changeStructureConductivity(O, f, er, ur, thk, sigma);
    end
    
    %% Public static function definitions (implemented in separate files)
    methods (Static, Access = public)
        [nodes, weights, errorWeights] = gaussKronrod(numSegs, a, b);
        [nodes, weights, errorWeights] = fejer2(orderN, a, b);
        [q] = integralVectorized(fun, a, b, options);
        
        [er, ur, thk] = validateStructure(f, er, ur, thk, options);
        [structureString, figureString] = printStructure(er, ur, thk, options);
    end

end

