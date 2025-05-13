classdef nLayerCircular < nLayerOpenEnded
    %Implementation of nLayerForward for circular waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by a circular waveguide looking into a multilayer structure. Note
    % that the units of all parameters should match that of the speed of
    % light specified by the "distanceUnitScale" and "frequencyUnitScale"
    % parameters (defaults are mm and GHz), both of which can be set in the
    % constructor.
    %
    % Example Usage:
    %   NL = nLayerCircular(0, maxN, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
    %   NL = nLayerCircular(m, maxN, waveguideR=3);
    %   NL = nLayerCircular(m, maxN, distanceUnitScale=1, waveguideR=3e-3);
    %   NL = nLayerCircular(m, maxN, prop1=val1, prop2=val2, ...);
    %
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, {}, thk);
    %   gam = NL.calculate(f, {}, ur, thk);
    %
    %
    % Author: Matt Dvorsky

    properties (Access=public, AbortSet)
        waveguideBand(1, 1) nLayer.circularWaveguideBand;   % Waveguide band.

        modeIndexM(1, 1) {mustBeInteger, mustBeNonnegative} = 0;    % Value of 'm' for considered TEmn or TMmn modes.
        maxModeIndexN(1, 1) {mustBeInteger, mustBePositive} = 1;    % Maximum value of 'n' for considered TEmn or TMmn modes.
    end
    properties (Dependent, Access=public, AbortSet)
        waveguideR(1, 1) {mustBePositive, mustBeFinite};    % Waveguide radius.
    end
    properties (Access=private)
        waveguideR_custom(1, 1);
    end

    %% Class Constructor
    methods
        function self = nLayerCircular(indexM, maxIndexN, classProperties)
            %Construct an instance of this class.
            % Inputs:
            %   indexM - Index 'm' of TEmn and TMmn modes to consider.
            %   maxIndexN - Highest index 'n' of TEmn and TMmn modes to consider.

            arguments
                indexM(1, 1) {mustBeInteger, mustBeNonnegative};
                maxIndexN(1, 1) {mustBeInteger, mustBePositive};
                classProperties.?nLayerCircular;
            end

            self.modeIndexM = indexM;
            self.maxModeIndexN = maxIndexN;

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                self.(propPairs{ii}) = propPairs{ii + 1};
            end

            switch self.modeSymmetryAxial
                case "TE"
                    self.frequencyRange = [3.8317, 7.0156] ./ (2*pi) ...
                        .* (self.speedOfLight ./ self.waveguideR);
                case "TM"
                    self.frequencyRange = [2.4048, 5.5201] ./ (2*pi) ...
                        .* (self.speedOfLight ./ self.waveguideR);
                case "None"
                    self.frequencyRange = [1.8412, 2.4048] ./ (2*pi) ...
                        .* (self.speedOfLight ./ self.waveguideR);
            end

            if isfield(classProperties, "frequencyRange")
                self.frequencyRange = classProperties.frequencyRange;
            end
        end
    end

    %% Class Functions
    methods (Access=protected)
        [waveguideModes] = defineWaveguideModes(self);
    end

    %% Class Setters
    methods
        function set.waveguideBand(self, newBand)
            self.waveguideBand = newBand;
            self.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.waveguideR(self, newR)
            self.waveguideR_custom = newR;
            self.waveguideBand = "Custom";
        end

        function set.modeIndexM(self, newMaxInd)
            self.modeIndexM = newMaxInd;
            self.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.maxModeIndexN(self, newMaxInd)
            self.maxModeIndexN = newMaxInd;
            self.shouldRegenerateWaveguideModeObjects = true;
        end
    end

    %% Class Getters
    methods
        function [r] = get.waveguideR(self)
            if self.waveguideBand == "Custom"
                r = self.waveguideR_custom;
                return;
            end
            [r] = self.waveguideBand.getDimensions(self.distanceUnitScale);
        end
    end

end

