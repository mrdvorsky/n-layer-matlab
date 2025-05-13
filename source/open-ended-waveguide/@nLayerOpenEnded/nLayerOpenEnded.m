classdef nLayerOpenEnded < nLayerForward
    %Implementation of nLayerForward for open-ended waveguides.
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
        waveguideModes(1, :) ...   % Array of "waveguideMode" objects, defining the properties of each mode.
            waveguideMode = waveguideMode.empty;
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
        [Smn] = calculate(self, f, er, ur, thk);
        [outputLabels] = getOutputLabels(self);
    end
    methods (Access=protected)
        [waveguideModes] = defineWaveguideModes(self, ...
            symmetryX, symmetryY, symmetryAxial);
    end
    methods (Access=private)
        [] = computeIntegralWeights(self, options);
        [krc, AhHat, AeHat] = computeAhat(self);
        [A] = computeA(self, f, er, ur, thk);
        [K] = computeK(self, f);
    end

    %% Class Constructor
    methods (Access=public)
        function self = nLayerOpenEnded(waveguideModes, classProperties)
            arguments
                waveguideModes(:, 1) waveguideMode = waveguideMode.empty;
                classProperties.?nLayerOpenEnded;
            end

            if ~isempty(waveguideModes)
                self.waveguideModes = waveguideModes;
            end

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                self.(propPairs{ii}) = propPairs{ii + 1};
            end
        end
    end

    %% Class Setters
    methods
        function set.frequencyRange(self, newFreqRange)
            self.frequencyRange_private = [min(newFreqRange), max(newFreqRange)];
            self.shouldRecomputeWeights = true;
        end

        function set.modeSymmetryX(self, newSym)
            self.modeSymmetryX_private = newSym;
            if strcmp(newSym, "PEC") && strcmp(self.modeSymmetryAxial, "TM")
                self.modeSymmetryAxial_private = "None";
            elseif strcmp(newSym, "PMC") && strcmp(self.modeSymmetryAxial, "TE")
                self.modeSymmetryAxial_private = "None";
            end
            self.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.modeSymmetryY(self, newSym)
            self.modeSymmetryY_private = newSym;
            if strcmp(newSym, "PEC") && strcmp(self.modeSymmetryAxial, "TM")
                self.modeSymmetryAxial_private = "None";
            elseif strcmp(newSym, "PMC") && strcmp(self.modeSymmetryAxial, "TE")
                self.modeSymmetryAxial_private = "None";
            end
            self.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.modeSymmetryAxial(self, newSym)
            self.modeSymmetryAxial_private = newSym;
            if strcmp(newSym, "TE")
                self.modeSymmetryX_private = "PEC";
                self.modeSymmetryY_private = "PEC";
            elseif strcmp(newSym, "TM")
                self.modeSymmetryX_private = "PMC";
                self.modeSymmetryY_private = "PMC";
            end
            self.shouldRegenerateWaveguideModeObjects = true;
        end
    end

    %% Class Getters
    methods
        function [waveguideModes] = get.waveguideModes(self)
            self.regenerateWaveguideModeObjects();
            waveguideModes = self.waveguideModes;
        end

        function [freq_range] = get.frequencyRange(self)
            freq_range = self.frequencyRange_private;
        end
        function [symX] = get.modeSymmetryX(self)
            symX = self.modeSymmetryX_private;
        end
        function [symY] = get.modeSymmetryY(self)
            symY = self.modeSymmetryY_private;
        end
        function [symAx] = get.modeSymmetryAxial(self)
            symAx = self.modeSymmetryAxial_private;
        end

        function [kc0] = get.mode_kc0(self)
            kc0 = [self.waveguideModes.kc0];
        end
        function [fc] = get.mode_fc0(self)
            fc = self.speedOfLight * self.mode_kc0 ./ (2*pi);
        end
        function [types] = get.modeTypes(self)
            types = [self.waveguideModes.modeType];
        end
        function [labels] = get.modeLabels(self)
            labels = [self.waveguideModes.modeLabel];
        end
        function [num] = get.numModes(self)
            num = numel(self.waveguideModes);
        end
        function [num] = get.numModes_TE(self)
            num = sum(strcmp(self.modeTypes, "TE"));
        end
        function [num] = get.numModes_TM(self)
            num = sum(strcmp(self.modeTypes, "TM"));
        end
        function [num] = get.numModes_Hybrid(self)
            num = sum(strcmp(self.modeTypes, "Hybrid"));
        end
    end

end

