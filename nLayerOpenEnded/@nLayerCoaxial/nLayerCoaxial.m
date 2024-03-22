classdef nLayerCoaxial < nLayerOpenEnded
    %NLAYERCOAXIAL Implementation of nLayerForward for coaxial waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by a rectangular waveguide looking into a multilayer structure. Note
    % that the units of all parameters should match that of the speed of
    % light specified by the speedOfLight parameter (default is mm GHz).
    %
    % Example Usage:
    %   NL = nLayerCoaxial(maxM, maxN, waveguideBand="N");
    %   NL = nLayerCoaxial(maxM, maxN, waveguideRi=0.44, waveguideRi=1);
    %   NL = nLayerCoaxial(maxM, maxN, distanceUnitScale=1, ...
    %       waveguideRi=0.44e-3, waveguideRi=1e-3);
    %   NL = nLayerCoaxial(maxM, maxN, waveguideBand="K", ...
    %       prop1=val1, prop2=val2, ...);
    %
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, [], thk);
    %   gam = NL.calculate(f, [], ur, thk);
    %   gam = NL.calculate(f, er, [], thk, BackingConductivity=sigma);
    %   gam = NL.calculate(f, er, [], thk, BackingConductivity=-inf);   % PMC 
    %
    % nLayerCoaxial Properties:
    %   waveguideRi - Coaxial inner radius.
    %   waveguideRo - Coaxial outer radius.
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
        waveguideBand(1, 1) coaxialConnectorType;       % Waveguide band.

        maxModeIndexM(1, 1) {mustBeInteger, mustBeNonnegative} = 1; % Maximum value of 'm' for considered TEmn or TMmn modes.
        maxModeIndexN(1, 1) {mustBeInteger, mustBeNonnegative} = 0; % Maximum value of 'n' for considered TEmn or TMmn modes.
    end
    properties (Dependent, Access=public, AbortSet)
        waveguideRi(1, 1) {mustBePositive, mustBeFinite};   % Coax inner radius.
        waveguideRo(1, 1) {mustBePositive, mustBeFinite};   % Coax outer radius.
    end
    properties (Access=private)
        waveguideRi_custom(1, 1);
        waveguideRo_custom(1, 1);
    end

    %% Class Functions
    methods (Access=protected)
        [modeStruct] = defineWaveguideModes(O);
    end

    %% Class Constructor
    methods
        function O = nLayerCoaxial(maxIndexM, maxIndexN, classProperties)
            %NLAYERCOAXIAL Construct an instance of this class.
            % Inputs:
            %   maxIndexM - Highest index 'm' of TEmn and TMmn modes to consider.
            %   maxIndexN - Highest index 'n' of TEmn and TMmn modes to consider.

            arguments
                maxIndexM(1, 1) {mustBeInteger, mustBeNonnegative};
                maxIndexN(1, 1) {mustBeInteger, mustBeNonnegative};
                classProperties.?nLayerCoaxial;
            end

            O.modeSymmetryAxial = "TM";

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                O.(propPairs{ii}) = propPairs{ii + 1};
            end

            O.maxModeIndexM = maxIndexM;
            O.maxModeIndexN = maxIndexN;

            rRatio = O.waveguideRi ./ O.waveguideRo;
            switch O.modeSymmetryAxial
                case "TM"
                    O.frequencyRange = [0, 0.99] ./ (2*pi) ...
                        .* (2.4048 * O.speedOfLight ./ O.waveguideRo ...
                        ./ (1 - rRatio));
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
        function set.waveguideRi(O, newA)
            O.waveguideRi_custom = newA;
            O.waveguideRo_custom = O.waveguideRo;
            O.waveguideBand = "Custom";
        end
        function set.waveguideRo(O, newB)
            O.waveguideRo_custom = newB;
            O.waveguideRi_custom = O.waveguideRi;
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
        function [ri] = get.waveguideRi(O)
            if O.waveguideBand == "Custom"
                ri = O.waveguideRi_custom;
                return;
            end
            [di, ~] = O.waveguideBand.getDimensions(O.distanceUnitScale);
            ri = 0.5 * di;
        end
        function [ro] = get.waveguideRo(O)
            if O.waveguideBand == "Custom"
                ro = O.waveguideRo_custom;
                return;
            end
            [~, do] = O.waveguideBand.getDimensions(O.distanceUnitScale);
            ro = 0.5 * do;
        end
    end

end

