classdef nLayerRectangular < nLayerOpenEnded
    %Implementation of nLayerForward for rectangular waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by a rectangular waveguide looking into a multilayer structure. Note
    % that the units of all parameters should match that of the speed of
    % light specified by the "distanceUnitScale" and "frequencyUnitScale"
    % parameters (defaults are mm and GHz), both of which can be set in the
    % constructor.
    %
    % Example Usage:
    %   NL = nLayerRectangular(maxM, maxN, waveguideBand="X");
    %   NL = nLayerRectangular(maxM, maxN, waveguideA=7.112, waveguideB=3.556);
    %   NL = nLayerRectangular(maxM, maxN, distanceUnitScale=1, ...
    %       waveguideA=7.112e-3, waveguideB=3.556e-3);
    %   NL = nLayerRectangular(maxM, maxN, waveguideBand="Ka", ...
    %       prop1=val1, prop2=val2, ...);
    %
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, {}, thk);
    %   gam = NL.calculate(f, {}, ur, thk);
    %
    %
    % Author: Matt Dvorsky

    properties (Access=public, AbortSet)
        waveguideBand(1, 1) nLayer.rectangularWaveguideBand;    % Waveguide band.

        maxModeIndexM(1, 1) {mustBeInteger, mustBeNonnegative} = 1; % Maximum value of 'm' for considered TEmn or TMmn modes.
        maxModeIndexN(1, 1) {mustBeInteger, mustBeNonnegative} = 0; % Maximum value of 'n' for considered TEmn or TMmn modes.
    end
    properties (Dependent, Access=public, AbortSet)
        waveguideA(1, 1) {mustBePositive, mustBeFinite};    % Waveguide broad dimension (along x-axis).
        waveguideB(1, 1) {mustBePositive, mustBeFinite};    % Waveguide narrow dimension (along y-axis).
    end
    properties (Access=private)
        waveguideA_custom(1, 1);
        waveguideB_custom(1, 1);
    end

    %% Class Constructor
    methods
        function self = nLayerRectangular(maxIndexM, maxIndexN, classProperties)
            %NLAYERRECTANGULAR Construct an instance of this class.
            % Inputs:
            %   maxIndexM - Highest index 'm' of TEmn and TMmn modes to consider.
            %   maxIndexN - Highest index 'n' of TEmn and TMmn modes to consider.

            arguments
                maxIndexM(1, 1) {mustBeInteger, mustBeNonnegative};
                maxIndexN(1, 1) {mustBeInteger, mustBeNonnegative};
                classProperties.?nLayerRectangular;
            end

            self.maxModeIndexM = maxIndexM;
            self.maxModeIndexN = maxIndexN;

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                self.(propPairs{ii}) = propPairs{ii + 1};
            end

            self.frequencyRange = [1.0, 2.0] ...
                .* (0.5 * self.speedOfLight ./ self.waveguideA);

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
        function set.waveguideA(self, newA)
            self.waveguideA_custom = newA;
            self.waveguideB_custom = self.waveguideB;
            self.waveguideBand = "Custom";
        end
        function set.waveguideB(self, newB)
            self.waveguideB_custom = newB;
            self.waveguideA_custom = self.waveguideA;
            self.waveguideBand = "Custom";
        end

        function set.maxModeIndexM(self, newMaxInd)
            self.maxModeIndexM = newMaxInd;
            self.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.maxModeIndexN(self, newMaxInd)
            self.maxModeIndexN = newMaxInd;
            self.shouldRegenerateWaveguideModeObjects = true;
        end
    end

    %% Class Getters
    methods
        function [a] = get.waveguideA(self)
            if self.waveguideBand == "Custom"
                a = self.waveguideA_custom;
                return;
            end
            [a, ~] = self.waveguideBand.getDimensions(self.distanceUnitScale);
        end
        function [b] = get.waveguideB(self)
            if self.waveguideBand == "Custom"
                b = self.waveguideB_custom;
                return;
            end
            [~, b] = self.waveguideBand.getDimensions(self.distanceUnitScale);
        end
    end

end

