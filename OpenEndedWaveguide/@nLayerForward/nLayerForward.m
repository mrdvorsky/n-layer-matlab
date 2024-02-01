classdef (Abstract) nLayerForward < matlab.mixin.Copyable & matlab.mixin.SetGetExactNames
    %NLAYERFORWARD Interface class for nLayer forward calculators.
    % This class serves as an interface definition for all nLayer forward
    % calculator objects. These objects take in a multilayer structure
    % definition and output the computed reflection coefficients (or
    % transmission coefficients).
    %
    % nLayerForward Properties:
    %   speedOfLight (299.792458) - Speed of light (default is mm GHz).
    %       Must match units of distance and frequency used.
    %   verbosity (0) - A value of 0 should suppress console output.
    %   checkStructureValues (true) - Flag used in the "verifyStructure"
    %       function. If true, this function will throw errors if
    %       non-physical values of er, ur, or thk are passed in.
    %
    % To use this class, subclass it. Any subclasses must, at a minimum,
    % implement the "calculate_impl" function and the "getOutputLabels"
    % function, whose functionalities are described below.
    % 
    % Abstract functions (should be implemented by subclass):
    %   calculate_impl - Should accept a vector of frequencies and the
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
    % Author: Matt Dvorsky
    
    properties (Access=public)
        speedOfLight(1, 1) {mustBePositive} = 299.792458;   % Speed of light (default is mm GHz).
        verbosity(1, 1) {mustBeNonnegative} = 0;            % Verbosity level. Set to zero for no console output.
        checkStructureValues(1, 1) logical = true;          % Whether to check ranges of er, ur, and thk.
    end
    
    %% Class Functions
    methods (Abstract, Access=protected)
        [gam] = calculate_impl(O, f, er, ur, thk, options);
    end
    methods (Abstract, Access=public)
        [outputLabels] = getOutputLabels(O);
    end
    methods (Access=public)
        [gam] = calculate(O, f, er, ur, thk);
        [er, ur, thk] = changeStructureConductivity(O, f, er, ur, thk, sigma);
    end
    methods (Static, Access=public)
        [] = checkClassProperties(O, classProperties);
    end

end

