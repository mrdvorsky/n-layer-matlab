classdef nLayerCircular < nLayerOpenEnded
    %NLAYERCIRCULAR Implementation of nLayerForward for circular waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by a rectangular waveguide looking into a multilayer structure. Note
    % that the units of all parameters should match that of the speed of
    % light specified by the speedOfLight parameter (default is mm GHz).
    %
    % Example Usage:
    %   NL = nLayerCircular(0, maxN, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
    %   NL = nLayerCircular(maxM, maxN, waveguideR=6);
    %   NL = nLayerCircular(maxM, maxN, distanceUnitScale=1, waveguideR=6e-3);
    %   NL = nLayerCircular(maxM, maxN, prop1=val1, prop2=val2, ...);
    %
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, [], thk);
    %   gam = NL.calculate(f, [], ur, thk);
    %   gam = NL.calculate(f, er, [], thk, BackingConductivity=sigma);
    %
    % nLayerRectangular Properties:
    %   waveguideR - Waveguide radius dimension.
    %   speedOfLight (299.792458) - Speed of light (default is mm GHz).
    %   verbosity (0) - Verbosity level. Set to 1 for basic command line
    %       output. Set to 2 for a per-frequency output.
    %   convergenceAbsTol (0.001) - Default tolerance for reflection
    %       coefficient calculations.
    %   checkStructureValues (true) - Flag used in the "verifyStructure"
    %       function. If true, this function will throw errors if
    %       non-physical values of er, ur, or thk are passed in.
    %
    % Author: Matt Dvorsky

    properties (Access=public, AbortSet)
        waveguideBand(1, 1) circularWaveguideBand;          % Waveguide band.

        maxModeIndexM(1, 1) {mustBeInteger, mustBeNonnegative} = 0; % Maximum value of 'm' for considered TEmn or TMmn modes.
        maxModeIndexN(1, 1) {mustBeInteger, mustBePositive} = 1;    % Maximum value of 'n' for considered TEmn or TMmn modes.
    end
    properties (Dependent, Access=public, AbortSet)
        waveguideR(1, 1) {mustBePositive, mustBeFinite};    % Waveguide radius.
    end
    properties (Access=private)
        waveguideR_custom(1, 1);
    end

    %% Class Functions
    methods (Access=protected)
        [modeStruct] = defineWaveguideModes(O);
    end

    %% Class Constructor
    methods
        function O = nLayerCircular(maxIndexM, maxIndexN, classProperties)
            %NLAYERCIRCULAR Construct an instance of this class.
            % Inputs:
            %   maxIndexM - Highest index 'm' of TEmn and TMmn modes to consider.
            %   maxIndexN - Highest index 'n' of TEmn and TMmn modes to consider.

            arguments
                maxIndexM(1, 1) {mustBeInteger, mustBeNonnegative};
                maxIndexN(1, 1) {mustBeInteger, mustBePositive};
            end
            arguments
                classProperties.?nLayerCircular;
            end

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                O.(propPairs{ii}) = propPairs{ii + 1};
            end

            O.maxModeIndexM = maxIndexM;
            O.maxModeIndexN = maxIndexN;

            switch O.modeSymmetryAxial
                case "TE"
                    O.frequencyRange = [1.01, 1.99] ./ (2*pi) ...
                        .* (3.8317 * O.speedOfLight ./ O.waveguideR);
                case "TM"
                    O.frequencyRange = [1.01, 1.99] ./ (2*pi) ...
                        .* (2.4048 * O.speedOfLight ./ O.waveguideR);
                case "None"
                    O.frequencyRange = [1.85, 2.39] ./ (2*pi) ...
                        .* (O.speedOfLight ./ O.waveguideR);
            end
            
            if O.numModes > 0
                O.excitationModeIndices = 1;
                O.receiveModeIndices = 1;
            end
        end
    end

    %% Class Setters
    methods
        function set.waveguideBand(O, newBand)
            O.waveguideBand = newBand;
            O.regenerateModeStructs();
        end
        function set.waveguideR(O, newR)
            O.waveguideR_custom = newR;
            O.waveguideBand = "Custom";
        end

        function set.maxModeIndexM(O, newMaxInd)
            O.maxModeIndexM = newMaxInd;
            O.regenerateModeStructs();
        end
        function set.maxModeIndexN(O, newMaxInd)
            O.maxModeIndexN = newMaxInd;
            O.regenerateModeStructs();
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

