classdef nLayerClamped < nLayerForward
    %NLAYERCLAMPED Implementation of nLayerForward for clamped waveguides.
    % This class can be used to calculate the reflection/transmission
    % coefficients of waveguides clamped on either side of an odd layered
    % or single layered structure.
    % Note that the units of all parameters should match that of
    % the speed of light specified by the speedOfLight parameter (default
    % is mm GHz).
    %
    % Example Usage:
    %   NL = nLayerCircular(0, 15, waveguideBand="X_TE01", modeSymmetryAxial="TE");
    %   NLC = nLayerClamped(NL);
    %
    %   S = NLC.calculate(f, er, ur, thk);
    %   S11 = S(:, 1, 1);
    %   S12 = S(:, 1, 2);
    %
    % nLayerClamped Properties:
    %   Uses same properties as as the passed nLayerForward object. 
    %
    % Author: Trent Moritz

    properties (GetAccess=public, SetAccess=public)
        outputIndices(:, 1) {mustBeInteger, mustBeInRange(outputIndices, 1, 4)} = [];   % List of enabled channel indices
        magCond(1, 1) {mustBePositive} = 10^12;   % Magnetic conductivity, used for PMC image
        NL(1,1); %nLayer object to be used for clamped setup.
    end

    %% Class Functions
    methods (Access=protected)
        [gam] = calculate_impl(O, f, er, ur, thk);
    end
    methods (Access=public)
        [outputLabels] = getOutputLabels(O);
    end

    %% Class constructor
    methods
        function O = nLayerClamped(NL, classProperties)
            %NLAYERCLAMPED Construct an instance of this class.
            % Example Usage:
            %   See example usage in main class documentation. Note that
            %   all public class properties can be specified as a named
            %   argument to the constructor (e.g., as "verbosity=1").
            
            arguments
                NL(1,1);
            end
            arguments (Repeating)
                classProperties;
            end

            %% Set Class Parameter Values
            O.NL = NL;

            if mod(numel(classProperties), 2) ~= 0
                error("Parameter and value arguments must come in pairs.");
            end

            if NL == 0
                error("Must pass nLayer object.")
            end

        end
    end

end

