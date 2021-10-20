classdef nLayerForward < handle
    %NLAYERFORWARD Interface class for nLayer forward calculators.
    % This class serves as an interface definition for all nLayer forward
    % calculator objects. These objects take in a multilayer structure
    % definition and output the computed reflection coefficients (or
    % transmission coefficients). It also contains several useful utility
    % functions.
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
    % "verbose" parameter value.
    %
    % Utility functions (implemented by this class):
    % gaussLegendre: Generates weights and nodes to perform Gaussian
    %   quadrature integration.
    % gaussKronrod: Generates weights and nodes to perform Gaussian
    %   quadrature integration.
    % integralVectorized: Routine to quickly perform integration of
    %   vectorized functions.
    % verifyStructure: Checks the validity of the multilayer
    %   structure and frequency definitions. This function should
    %   be called at the beginning of calculate.
    
    properties (GetAccess = public, SetAccess = public)       
        thkMin = inf;               % Minimum layer thickness.
        thkMax = 0;                 % Maximum layer thickness.
        erMin = 1;                  % Minimum layer complex permittivity.
        erMax = inf - 1j*inf;       % Maximum layer complex permittivity.
        urMin = 1;                  % Minimum layer complex permeability.
        urMax = inf - 1j*inf;       % Maximum layer complex permeability.
        
        c = 299.792458;             % Speed of light (mm GHz). Defines the units used for all calculations.
        
        verbosity = 0;              % Verbosity Level. Set to 1 for command line output.
    end
    
    %% Public member function definitions (virtual, not implemented)
    methods (Access = public)
        gam = calculate(O, f, er, ur, thk);
    end
    
    %% Public member function definitions
    methods (Access = public)
        function [f, er, ur, thk] = verifyStructure(O, f, er, ur, thk)
            
            % Check dimensions
        end
    end
    
    %% Public static function definitions (implemented in separate files)
    methods (Static, Access = public)
        [nodes, weights, errorWeights] = gaussKronrod(numSegs, a, b);
        [nodes, weights] = gaussLegendre(orderN, a, b);
        [q] = integralVectorized(fun, a, b, options);
        [q, nodes_out, weights_out] = integralWeightsAndNodes(fun, a, b, options);
    end

end

