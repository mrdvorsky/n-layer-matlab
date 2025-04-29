classdef nLayerOpenEnded < nLayerForward
    %NLAYEROPENENDED Implementation of nLayerForward for open-ended waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by an arbitrary waveguide(s) looking into a multilayer structure.
    % This class is meant to either be subclassed, where the subclass
    % defines the modes to be considered, or an array of
    % "nLayer.waveguideMode" objects can be passed into the constructor.
    %
    % Note: The units of all parameters should match that of the speed of
    % light specified by the "speedOfLight", "distanceUnitScale", and
    % "frequencyUnitScale" parameters (defaults are mm and GHz).
    %
    % Author: Matt Dvorsky

    %% Public Properties
    properties (GetAccess=public, SetAccess=protected)
        waveguideModes(1, :) ...   % Array of "nLayer.waveguideMode" objects, defining the properties of each mode.
            nLayer.waveguideMode = nLayer.waveguideMode.empty;
    end
    properties (Dependent, Access=public, AbortSet)
        frequencyRange(1, :) {mustBeNonnegative, mustBeFinite}; % Operating frequency range of the object.
        
        modeSymmetryX(1, 1) string ...          % Symmetry of reflection about the x-axis.
            {mustBeMember(modeSymmetryX, ["PEC", "PMC", "None"])};
        modeSymmetryY(1, 1) string ...          % Symmetry of reflection about the y-axis.
            {mustBeMember(modeSymmetryY, ["PEC", "PMC", "None"])};
        modeSymmetryAxial(1, 1) string ...      % Axial symmetry.
            {mustBeMember(modeSymmetryAxial, ["TE", "TM", "None"])};
    end
    properties (Access=public)
        waveguideEr(1, 1) {nLayer.mustBeErUrCallable} = 1;  % Array of filled waveguide permittivity values or function handle.
        waveguideUr(1, 1) {nLayer.mustBeErUrCallable} = 1;  % Array of filled waveguide permeability values or function handle.

        excitationModeIndices(1, :) {mustBeInteger, mustBePositive} = [1];  % Array of indices of excitation modes (i.e., 'n' in Smn).
        receiveModeIndices(1, :) {mustBeInteger, mustBePositive} = [1];     % Array of indices of receiving modes (i.e., 'm' in Smn).
    end
    properties(Dependent, GetAccess=public)
        mode_kc0;           % Array of free-space cutoff wavenumbers (kc0) for each mode.
        mode_fc0;           % Array of free-space cutoff frequencies for each mode.
        modeTypes;          % String array of mode type for each mode ("TE", "TM", or "Hybrid").
        modeLabels;         % String array of labels for each mode.

        numModes;           % Total number of modes considered (TE + TM + Hybrid).
        numModes_TE;        % Number of TE modes considered.
        numModes_TM;        % Number of TM modes considered.
        numModes_Hybrid;    % Number of Hybrid modes considered.
    end

    %% Private Properties
    properties (Hidden, Access=protected)
        frequencyRange_private(1, 2) {mustBeNonnegative, mustBeFinite} = [1, 2];
        modeSymmetryX_private(1, 1) string ...
            {mustBeMember(modeSymmetryX_private, ["PEC", "PMC", "None"])} = "PEC";
        modeSymmetryY_private(1, 1) string ...
            {mustBeMember(modeSymmetryY_private, ["PEC", "PMC", "None"])} = "PMC";
        modeSymmetryAxial_private(1, 1) string ...
            {mustBeMember(modeSymmetryAxial_private, ["TE", "TM", "None"])} = "None";
    end
    properties (Hidden, Access=protected)
        fixed_kr;       % Fixed-point integral coordindates kr.
        fixed_Ah;       % Fixed-point integral weights for Ah.
        fixed_Ae;       % Fixed-point integral weights for Ae.
        
        shouldRecomputeWeights(1, 1) logical = true;    % Flag to recompute integral weights.
        shouldRegenerateWaveguideModeObjects(1, 1) logical = true;  % Flag to regenerate waveguideMode objects.
    end
    properties (Access=public)
        integral_pointsKrc = {50, 50, 50, 50, 50};
        integral_pointsKr = {1024, 4096, 4096, 4096, 4096};
        integral_pointsPhi = 64;
    end

    %% Class Functions
    methods (Access=public)
        [Smn] = calculate(O, f, er, ur, thk);
        [outputLabels] = getOutputLabels(O);
    end
    methods (Access=protected)
        [waveguideModes] = defineWaveguideModes(O, ...
            symmetryX, symmetryY, symmetryAxial);
    end
    methods (Access=private)
        [] = computeIntegralWeights(O, options);
        [krc, AhHat, AeHat] = computeAhat(O);
        [A] = computeA(O, f, er, ur, thk);
        [K] = computeK(O, f);
    end

    %% Class Constructor
    methods (Access=public)
        function O = nLayerOpenEnded(waveguideModes, classProperties)
            arguments
                waveguideModes(:, 1) nLayer.waveguideMode = nLayer.waveguideMode.empty;
                classProperties.?nLayerOpenEnded;
            end

            if ~isempty(waveguideModes)
                O.waveguideModes = waveguideModes;
            end

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                O.(propPairs{ii}) = propPairs{ii + 1};
            end
        end
    end

    %% Class Setters
    methods
        function set.frequencyRange(O, newFreqRange)
            O.frequencyRange_private = [min(newFreqRange), max(newFreqRange)];
            O.shouldRecomputeWeights = true;
        end

        function set.modeSymmetryX(O, newSym)
            O.modeSymmetryX_private = newSym;
            if strcmp(newSym, "PEC") && strcmp(O.modeSymmetryAxial, "TM")
                O.modeSymmetryAxial_private = "None";
            elseif strcmp(newSym, "PMC") && strcmp(O.modeSymmetryAxial, "TE")
                O.modeSymmetryAxial_private = "None";
            end
            O.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.modeSymmetryY(O, newSym)
            O.modeSymmetryY_private = newSym;
            if strcmp(newSym, "PEC") && strcmp(O.modeSymmetryAxial, "TM")
                O.modeSymmetryAxial_private = "None";
            elseif strcmp(newSym, "PMC") && strcmp(O.modeSymmetryAxial, "TE")
                O.modeSymmetryAxial_private = "None";
            end
            O.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.modeSymmetryAxial(O, newSym)
            O.modeSymmetryAxial_private = newSym;
            if strcmp(newSym, "TE")
                O.modeSymmetryX_private = "PEC";
                O.modeSymmetryY_private = "PEC";
            elseif strcmp(newSym, "TM")
                O.modeSymmetryX_private = "PMC";
                O.modeSymmetryY_private = "PMC";
            end
            O.shouldRegenerateWaveguideModeObjects = true;
        end
    end

    %% Class Getters
    methods
        function [waveguideModes] = get.waveguideModes(O)
            O.regenerateWaveguideModeObjects();
            waveguideModes = O.waveguideModes;
        end

        function [freq_range] = get.frequencyRange(O)
            freq_range = O.frequencyRange_private;
        end
        function [symX] = get.modeSymmetryX(O)
            symX = O.modeSymmetryX_private;
        end
        function [symY] = get.modeSymmetryY(O)
            symY = O.modeSymmetryY_private;
        end
        function [symAx] = get.modeSymmetryAxial(O)
            symAx = O.modeSymmetryAxial_private;
        end

        function [kc0] = get.mode_kc0(O)
            kc0 = [O.waveguideModes.kc0];
        end
        function [fc] = get.mode_fc0(O)
            fc = O.speedOfLight * O.mode_kc0 ./ (2*pi);
        end
        function [types] = get.modeTypes(O)
            types = [O.waveguideModes.modeType];
        end
        function [labels] = get.modeLabels(O)
            labels = [O.waveguideModes.modeLabel];
        end
        function [num] = get.numModes(O)
            num = numel(O.waveguideModes);
        end
        function [num] = get.numModes_TE(O)
            num = sum(strcmp(O.modeTypes, "TE"));
        end
        function [num] = get.numModes_TM(O)
            num = sum(strcmp(O.modeTypes, "TM"));
        end
        function [num] = get.numModes_Hybrid(O)
            num = sum(strcmp(O.modeTypes, "Hybrid"));
        end
    end

end

