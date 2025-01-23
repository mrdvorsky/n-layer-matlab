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

    properties (GetAccess=public, SetAccess=protected)
        waveguideModes(1, :) ...   % Array of "nLayer.waveguideMode" objects, defining the properties of each mode.
            nLayer.waveguideMode = nLayer.waveguideMode.empty;
    end
    properties (Access=public, AbortSet)
        frequencyRange(1, :) {mustBeNonnegative, mustBeFinite} = [1, 2];   % Operating frequency range of the object.
        
        modeSymmetryX string ...            % Symmetry of reflection about the x-axis.
            {mustBeMember(modeSymmetryX, ["PEC", "PMC", "None"])} = "PEC";
        modeSymmetryY string ...            % Symmetry of reflection about the y-axis.
            {mustBeMember(modeSymmetryY, ["PEC", "PMC", "None"])} = "PMC";
        modeSymmetryAxial string ...        % Axial symmetry.
            {mustBeMember(modeSymmetryAxial, ["TE", "TM", "None"])} = "None";

        waveguideEr(1, 1) {nLayer.mustBeErUrCallable} = 1;  % Array of filled waveguide permittivity values for each mode.
        waveguideUr(1, 1) {nLayer.mustBeErUrCallable} = 1;  % Array of filled waveguide permeability values for each mode.

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
        [gam] = calculate(O, f, er, ur, thk);
        [outputLabels] = getOutputLabels(O);
    end
    methods (Access=protected)
        [modeStructs] = defineWaveguideModes(O, ...
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
        function O = nLayerOpenEnded(modeStructs)
            arguments
                modeStructs(:, 1) nLayer.waveguideMode = nLayer.waveguideMode.empty;
            end

            if ~isempty(modeStructs)
                O.waveguideModes = modeStructs;
            end
        end
    end

    %% Class Setters
    methods
        function set.frequencyRange(O, newFreqRange)
            O.frequencyRange = [min(newFreqRange), max(newFreqRange)];
            O.shouldRecomputeWeights = true;        %#ok<MCSUP>
        end

        function set.modeSymmetryX(O, newSym)
            O.modeSymmetryX = newSym;
            O.shouldRegenerateWaveguideModeObjects = true;  %#ok<MCSUP>
        end
        function set.modeSymmetryY(O, newSym)
            O.modeSymmetryY = newSym;
            O.shouldRegenerateWaveguideModeObjects = true;  %#ok<MCSUP>
        end
        function set.modeSymmetryAxial(O, newSym)
            O.modeSymmetryAxial = newSym;
            if strcmp(newSym, "TE")
                O.modeSymmetryX = "PEC";        %#ok<MCSUP>
                O.modeSymmetryY = "PEC";        %#ok<MCSUP>
            elseif strcmp(newSym, "TM")
                O.modeSymmetryX = "PMC";        %#ok<MCSUP>
                O.modeSymmetryY = "PMC";        %#ok<MCSUP>
            end
            O.shouldRegenerateWaveguideModeObjects = true;  %#ok<MCSUP>
        end
    end

    %% Class Getters
    methods
        function [waveguideModes] = get.waveguideModes(O)
            O.regenerateWaveguideModeObjects();
            waveguideModes = O.waveguideModes;
        end

        function [kc0] = get.mode_kc0(O)
            kc0 = [O.waveguideModes.CutoffWavenumber];
        end
        function [fc] = get.mode_fc0(O)
            fc = O.speedOfLight * O.mode_kc0 ./ (2*pi);
        end
        function [types] = get.modeTypes(O)
            types = [O.waveguideModes.ModeType];
        end
        function [labels] = get.modeLabels(O)
            labels = [O.waveguideModes.ModeLabel];
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

