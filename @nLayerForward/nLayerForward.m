classdef (Abstract) nLayerForward < handle
    %NLAYERFORWARD Interface class for nLayer forward calculators.
    % This class serves as an interface definition for all nLayer forward
    % calculator objects. These objects take in a multilayer structure
    % definition and output the computed reflection coefficients (or
    % transmission coefficients). It also contains several useful utility
    % functions.
    %
    % nLayerForward Properties:
    %   c (299.792458) - Speed of light (mm GHz). Must match units of
    %       distance and frequency used.
    %   verbosity (0) - A value of 0 should suppress console output.
    %   checkStructureValues (true) - Flag used in the "verifyStructure"
    %       function. If true, this function will throw errors if
    %       non-physical values of er, ur, or thk are passed in.
    %
    % To use this class, subclass it. Any subclasses must, at a minimum,
    % implement the "calculate" function, which should take a vector of
    % frequencies and the multilayer structure definition at each
    % frequency, and should output the calculated reflection (and/or
    % transmission) coefficient(s) at each frequency.
    %
    % The implementation of "calculate" should utilize the parameter "c"
    % for defining units (default mm GHz), the "verifyStructure" function
    % for checking the multilayer structure definition, and observe the
    % "verbosity" parameter value.
    %
    % Utility functions (implemented by this class):
    %   fejer2: Generates weights and nodes to perform Fejer Type II
    %       quadrature integration.
    %   gaussKronrod: Generates weights and nodes to perform Gaussian
    %       quadrature integration.
    %   integralVectorized: Routine to quickly perform adaptive integration
    %       of vectorized functions.
    %   integralWeightsAndNodes: Routine to perform integration of
    %       functions, and return the final weights and nodes used.
    %   verifyStructure: Checks the validity of the multilayer structure
    %       and frequency definitions. This function should be called at
    %       the beginning of calculate.
    %
    % Author: Matt Dvorsky
    
    properties (GetAccess = public, SetAccess = public)       
        c = 299.792458;
        verbosity = 0;
        checkStructureValues = true;
    end
    
    %% Virtual Public member function definitions
    methods (Abstract, Access = public)
        gam = calculate(O, f, er, ur, thk);
    end
    
    %% Public member function definitions (implemented in separate files)
    methods (Access = public)
        gam = calculateWithConductivity(O, f, er, ur, thk, sigma, varargin);
    end
    
    %% Public static function definitions (implemented in separate files)
    methods (Static, Access = public)
        [nodes, weights, errorWeights] = gaussKronrod(numSegs, a, b);
        [nodes, weights, errorWeights] = fejer2(orderN, a, b);
        [q] = integralVectorized(fun, a, b, options);
        [q, nodes_out, weights_out] = integralWeightsAndNodes(fun, a, b, options);
        
        [f, er, ur, thk] = verifyStructure(f, er, ur, thk, options);
        [er, ur, thk] = changeStructureConductivity(f, er, ur, thk, sigma);
    end

end

