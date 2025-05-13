classdef (Abstract) nLayerForward < matlab.mixin.Copyable
    %Interface class for nLayer forward calculators.
    % This class serves as an interface definition for all nLayer forward
    % calculator objects. These objects take in a multilayer structure
    % definition and output the computed reflection coefficients (or
    % transmission coefficients).
    %
    % To use this class, subclass it. Any subclasses must, at a minimum,
    % implement the "calculate" function and the "getOutputLabels"
    % function, whose functionalities are described below.
    %
    % Abstract functions (should be implemented by subclass):
    %   calculate - Should accept a vector of frequencies and the
    %       multilayer structure definition at each frequency, and should
    %       output the calculated reflection (and/or transmission)
    %       coefficient(s) at each frequency. The first column of the
    %       output should have the same length as the number of
    %       frequencies. Additional dimensions can be anything.
    %   getOutputLabels - Should return a vector of strings labeling each
    %       output channel. If "calculateGamma" returns gam, then the nth
    %       element of this output vector should describe gam(:, n).
    %
    % The implementation of "calculate" should utilize the parameters
    % "speedOfLight", "distanceUnitScale", and "frequencyUnitScale" for
    % defining units (default mm GHz), and should observe the "verbosity"
    % parameter value.
    %
    % Author: Matt Dvorsky

    properties (Access=public)
        verbosity(1, 1) {mustBeNonnegative} = 0;            % Verbosity level. Set to zero for no console output.
        checkStructureValues(1, 1) logical = true;          % Whether to check ranges of er, ur, and thk.
    end
    properties (GetAccess=public, SetAccess=immutable)
        distanceUnitScale(1, 1) {mustBePositive} = 0.001;   % All input distances will be scaled by this (default is for mm).
        frequencyUnitScale(1, 1) {mustBePositive} = 1e9;    % All input frequencies will be scaled by this (default is for GHz).
    end
    properties (Dependent, GetAccess=public)
        speedOfLight(1, 1);     % Speed of light, with units based on "distanceUnitScale" and "frequencyUnitScale".
    end

    %% Class Getters
    methods
        function [speedOfLight] = get.speedOfLight(self)
            speedOfLight = 299792458 ...
                ./ (self.distanceUnitScale * self.frequencyUnitScale);
        end
    end

    %% Class Functions
    methods (Abstract, Access=public)
        [gam] = calculate(self, f, er, ur, thk, options);
        [outputLabels] = getOutputLabels(self);
    end

end

