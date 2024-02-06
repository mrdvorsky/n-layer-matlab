classdef nLayerRectangular < nLayerOpenEnded
    %NLAYERRECTANGULAR Implementation of nLayerForward for rectangular waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by a rectangular waveguide looking into a multilayer structure. Note
    % that the units of all parameters should match that of the speed of
    % light specified by the speedOfLight parameter (default is mm GHz).
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
    %   gam = NL.calculate(f, er, [], thk);
    %   gam = NL.calculate(f, [], ur, thk);
    %   gam = NL.calculate(f, er, [], thk, BackingConductivity=sigma);
    %
    % nLayerRectangular Properties:
    %   waveguideA - Waveguide broad dimension.
    %   waveguideB - Waveguide narrow dimension.
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
        waveguideBand(1, 1) rectangularWaveguideBand;       % Waveguide band.

        maxModeIndexM(1, 1) {mustBeInteger, mustBeNonnegative} = 1; % Maximum value of 'm' for considered TEmn or TMmn modes.
        maxModeIndexN(1, 1) {mustBeInteger, mustBeNonnegative} = 0; % Maximum value of 'n' for considered TEmn or TMmn modes.
    end
    properties (Dependent, Access=public, AbortSet)
        waveguideA(1, 1) {mustBePositive, mustBeFinite};    % Waveguide broad dimension.
        waveguideB(1, 1) {mustBePositive, mustBeFinite};    % Waveguide narrow dimension.
    end
    properties (Access=private)
        waveguideA_custom(1, 1);
        waveguideB_custom(1, 1);
    end

    %% Class Functions
    methods (Access=protected)
        [modeStruct] = defineWaveguideModes(O);
    end

    %% Class Constructor
    methods
        function O = nLayerRectangular(maxIndexM, maxIndexN, classProperties)
            %NLAYERRECTANGULAR Construct an instance of this class.
            % Inputs:
            %   maxIndexM - Highest index 'm' of TEmn and TMmn modes to consider.
            %   maxIndexN - Highest index 'n' of TEmn and TMmn modes to consider.

            arguments
                maxIndexM(1, 1) {mustBeInteger, mustBePositive};
                maxIndexN(1, 1) {mustBeInteger, mustBeNonnegative};
            end
            arguments
                classProperties.?nLayerRectangular;
            end

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                O.(propPairs{ii}) = propPairs{ii + 1};
            end

            O.maxModeIndexM = maxIndexM;
            O.maxModeIndexN = maxIndexN;
        end
    end

    %% Class Setters
    methods
        function set.waveguideBand(O, newBand)
            O.waveguideBand = newBand;
            O.regenerateModeStructs();
        end
        function set.waveguideA(O, newA)
            O.waveguideA_custom = newA;
            O.waveguideB_custom = O.waveguideB;
            O.waveguideBand = "Custom";
        end
        function set.waveguideB(O, newB)
            O.waveguideB_custom = newB;
            O.waveguideA_custom = O.waveguideA;
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
        function a = get.waveguideA(O)
            if O.waveguideBand == "Custom"
                a = O.waveguideA_custom;
                return;
            end
            [a, ~] = O.waveguideBand.getDimensions(O.distanceUnitScale);
        end
        function b = get.waveguideB(O)
            if O.waveguideBand == "Custom"
                b = O.waveguideB_custom;
                return;
            end
            [~, b] = O.waveguideBand.getDimensions(O.distanceUnitScale);
        end
    end

end

