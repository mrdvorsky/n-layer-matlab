classdef nLayerFilledRectangular_old < nLayerForward
    %Implementation of nLayerForward for filled rectangular waveguides.
    % This class can be used to calculate the reflection/transmission
    % coefficients of a rectangular waveguide filled with a multilayer
    % structure. Note that the units of all parameters should match that of
    % the speed of light specified by the speedOfLight parameter (default
    % is mm GHz).
    %
    % Example Usage:
    %   NL = nLayerFilledRectangular(TE_m, TE_n, waveguideBand="x");
    %
    %   S = NL.calculate(f, er, ur, thk);
    %   S11 = S(:, 1, 1);
    %   S12 = S(:, 1, 2);
    %
    %   S = NL.calculate(f, er, [], thk, BackingConductivity=sigma);
    %   S11 = S(:, 1, 1);
    %
    % nLayerFilledRectangular Properties:
    %   speedOfLight (299.792458) - Speed of light (default is mm GHz).
    %   verbosity (0) - Verbosity level. Set to 1 for basic command line
    %       output. Set to 2 for a per-frequency output.
    %   checkStructureValues (true) - Flag used in the "verifyStructure"
    %       function. If true, this function will throw errors if
    %       non-physical values of er, ur, or thk are passed in.
    %
    % Author: Trent Moritz

    properties (GetAccess=public, SetAccess=public)
        outputIndices(:, 1) {mustBeInteger, mustBeInRange(outputIndices, 1, 4)} = [];   % List of enabled channel indices
        waveguideA(1, 1) {mustBePositive} = 1;                  % Waveguide broad dimension.
        waveguideB(1, 1) {mustBePositive} = 0.5;                % Waveguide narrow dimension.
        modeTE_m(1,1) {mustBeInteger, mustBeNonnegative} = 1;   %TE mode index m
        modeTE_n(1,1) {mustBeInteger, mustBeNonnegative} = 0;   %TE mode index n
    end
    properties (GetAccess=public, SetAccess=private)
        waveguideBand = "";     % Waveguide band identifier.
    end

    %% Class Functions
    methods (Access=public)
        [gam] = calculate(O, f, er, ur, thk);
        [outputLabels] = getOutputLabels(O);
        [waveguideA, waveguideB] = setWaveguideBand(O, band, options);
        setMode(O, modeTE_m, modeTE_n);
    end

    %% Class constructor
    methods
        function O = nLayerFilledRectangular_old(modeTE_m, modeTE_n, classProperties)
            %NLAYERFILLEDRECTANGULAR Construct an instance of this class.
            % Example Usage:
            %   See example usage in main class documentation. Note that
            %   all public class properties can be specified as a named
            %   argument to the constructor (e.g., as "verbosity=1").
            
            arguments
                modeTE_m(1,1) {mustBeInteger, mustBeNonnegative} = 1;
                modeTE_n(1,1) {mustBeInteger, mustBeNonnegative} = 0;
            end
            arguments (Repeating)
                classProperties;
            end

            %% Set Class Parameter Values
            if mod(numel(classProperties), 2) ~= 0
                error("Parameter and value arguments must come in pairs.");
            end
            for ii = 1:2:numel(classProperties)
                O.(classProperties{ii}) = classProperties{ii + 1};
            end
            if isempty(O.modeTE_m)
                O.setMode(modeTE_m, modeTE_n);
            end
            if strlength(O.waveguideBand) > 0
                O.setWaveguideBand(O.waveguideBand);
            end
            O.setMode(modeTE_m,modeTE_n);
        end
    end

end

