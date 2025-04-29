classdef nLayerPlaneWave < nLayerForward
    %NLAYERPLANEWAVE Implementation of nLayerForward for plane waves.
    % This class can be used to calculate the reflection/transmission
    % coefficients of a plane-wave incident upon a flat multi-layered
    % structure. Note that the units of all parameters should match that
    % of the speed of light specified by the speedOfLight parameter
    % (default is mm GHz).
    %
    % Example Usage:
    %   NL = nLayerPlaneWave("TE");
    %
    %   S = NL.calculate(f, er, ur, thk);
    %   S11 = S(:, 1, 1);
    %   S12 = S(:, 1, 2);
    %
    % Author: Matt Dvorsky

    properties (Access=protected)
        modeType(1, 1) string {mustBeMember(modeType, ["TE", "TM"])};
    end

    %% Class Functions
    methods (Access=public)
        [gam] = calculate(O, f, er, ur, thk);
        [outputLabels] = getOutputLabels(O);
    end

    %% Class constructor
    methods
        function O = nLayerPlaneWave(modeType, classProperties)
            %NLAYERPLANEWAVE Construct an instance of this class.
            
            arguments
                modeType(1, 1) string {mustBeMember(modeType, ["TE", "TM"])}
            end
            arguments (Repeating)
                classProperties.?nLayerPlaneWave;
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

