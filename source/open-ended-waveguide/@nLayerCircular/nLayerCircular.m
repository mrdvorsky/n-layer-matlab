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
        function O = nLayerCircular(indexM, maxIndexN, classProperties)
            %Construct an instance of this class.
            % Inputs:
            %   indexM - Index 'm' of TEmn and TMmn modes to consider.
            %   maxIndexN - Highest index 'n' of TEmn and TMmn modes to consider.

            arguments
                indexM(1, 1) {mustBeInteger, mustBeNonnegative};
                maxIndexN(1, 1) {mustBeInteger, mustBePositive};
                classProperties.?nLayerCircular;
            end

            O.modeIndexM = indexM;
            O.maxModeIndexN = maxIndexN;

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                O.(propPairs{ii}) = propPairs{ii + 1};
            end

            switch O.modeSymmetryAxial
                case "TE"
                    O.frequencyRange = [3.8317, 7.0156] ./ (2*pi) ...
                        .* (O.speedOfLight ./ O.waveguideR);
                case "TM"
                    O.frequencyRange = [2.4048, 5.5201] ./ (2*pi) ...
                        .* (O.speedOfLight ./ O.waveguideR);
                case "None"
                    O.frequencyRange = [1.8412, 2.4048] ./ (2*pi) ...
                        .* (O.speedOfLight ./ O.waveguideR);
            end

            if isfield(classProperties, "frequencyRange")
                O.frequencyRange = classProperties.frequencyRange;
            end
        end
    end

    %% Class Functions
    methods (Access=protected)
        [waveguideModes] = defineWaveguideModes(O);
    end

    %% Class Setters
    methods
        function set.waveguideBand(O, newBand)
            O.waveguideBand = newBand;
            O.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.waveguideR(O, newR)
            O.waveguideR_custom = newR;
            O.waveguideBand = "Custom";
        end

        function set.modeIndexM(O, newMaxInd)
            O.modeIndexM = newMaxInd;
            O.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.maxModeIndexN(O, newMaxInd)
            O.maxModeIndexN = newMaxInd;
            O.shouldRegenerateWaveguideModeObjects = true;
        end
    end

    %% Class Getters
    methods
        function [r] = get.waveguideR(O)
            if O.waveguideBand == "Custom"
                r = O.waveguideR_custom;
                return;
            end
            [r] = O.waveguideBand.getDimensions(O.distanceUnitScale);
        end
    end

end

