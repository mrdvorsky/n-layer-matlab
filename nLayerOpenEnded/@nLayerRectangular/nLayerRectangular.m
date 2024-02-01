classdef nLayerRectangular < nLayerOpenEnded
    %NLAYERRECTANGULAR Implementation of nLayerForward for rectangular waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by a rectangular waveguide looking into a multilayer structure. Note
    % that the units of all parameters should match that of the speed of
    % light specified by the speedOfLight parameter (default is mm GHz).
    %
    % Example Usage:
    %   NL = nLayerRectangular(maxM, maxN, waveguideBand="x");
    %   NL = nLayerRectangular(maxM, maxN, waveguideA=7.112, waveguideB=3.556);
    %   NL = nLayerRectangular(maxM, maxN, speedOfLight=299.79e6, ...
    %       waveguideA=7.112e-3, waveguideB=3.556e-3);
    %   NL = nLayerRectangular(maxM, maxN, waveguideBand="ka", ...
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
    %   modes_TE - List of TEmn modes to consider (in rows of [m, n]).
    %   modes_TM - List of TMmn modes to consider (in rows of [m, n]).
    %   verbosity (0) - Verbosity level. Set to 1 for basic command line
    %       output. Set to 2 for a per-frequency output.
    %   convergenceAbsTol (0.001) - Default tolerance for reflection
    %       coefficient calculations.
    %   checkStructureValues (true) - Flag used in the "verifyStructure"
    %       function. If true, this function will throw errors if
    %       non-physical values of er, ur, or thk are passed in.
    %
    % Author: Matt Dvorsky

    properties (GetAccess=public, SetAccess=public)
        waveguideA(1, 1) {mustBePositive} = 1;              % Waveguide broad dimension.
        waveguideB(1, 1) {mustBePositive} = 0.5;            % Waveguide narrow dimension.
        modes_TE(:, 2) {mustBeInteger, mustBeNonnegative};  % List of TEmn modes in rows of [m, n].
        modes_TM(:, 2) {mustBeInteger, mustBeNonnegative};  % List of TMmn modes in rows of [m, n].
    end
    properties (GetAccess=public, SetAccess=private)
        waveguideBand = "";         % Waveguide band identifier.
        modeSymmetryX = "Even";     % Mode symmetry specification along x-axis.
        modeSymmetryY = "Odd";      % Mode symmetry specification along y-axis.
    end

    %% Class Functions
    methods (Access=protected)
        [modeStruct] = defineWaveguideModes(O);
    end
    methods (Access=public)
        [modesTE, modesTM] = setModes(O, maxM, maxN);
        [waveguideA, waveguideB] = setWaveguideBand(O, band, options);
    end

    %% Class Constructor
    methods
        function O = nLayerRectangular(maxM, maxN, classProperties)
            %NLAYERRECTANGULAR Construct an instance of this class.
            % Inputs:
            %   maxM - Highest index m of TEmn and TMmn modes to consider.
            %   maxN - Highest index n of TEmn and TMmn modes to consider.

            arguments
                maxM(1, 1) {mustBeInteger, mustBePositive};
                maxN(1, 1) {mustBeInteger, mustBeNonnegative};
            end
            arguments (Repeating)
                classProperties;
            end

            %% Set Class Parameter Values
            O.checkClassProperties(O, classProperties);
            for ii = 1:2:numel(classProperties)
                set(O, classProperties{ii}, classProperties{ii + 1});
            end

            if isempty(O.modes_TE)
                O.setModes(maxM, maxN);
            end

            if strlength(O.waveguideBand) > 0
                O.setWaveguideBand(O.waveguideBand);
            end

            O.recomputeInterpolants();
        end
    end

end

