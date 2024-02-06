classdef nLayerOpenEnded < nLayerForward
    %NLAYEROPENENDED Implementation of nLayerForward for open-ended waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by an arbitrary waveguide(s) looking into a multilayer structure.
    % This class is meant to either be subclassed, where the subclass
    % defines the modes to be considered, or an array of modeStructs can be
    % passed into the constructor. Alternatively, an array of
    % nLayerOpenEnded objects can be passed into the constructor, which
    % will results in a new object considered all modes to be created.
    %
    % Note: The units of all parameters should match that of the speed of
    % light specified by the speedOfLight parameter (default is mm GHz).
    %
    % nLayerOpenEnded Properties:
    %   speedOfLight (299.792458) - Speed of light (default is mm GHz).
    %   verbosity (0) - Verbosity level. Set to 1 for basic command line
    %       output. Set to 2 for a per-frequency output.
    %
    % Author: Matt Dvorsky

    properties (GetAccess=public, SetAccess=public)
        modeStructs(1, :);      % Array of waveguideMode objects, defining the properties of each mode.
        frequencyRange(1, :) {mustBePositive, mustBeFinite} = [1, 2];   % Frequency range that the object will be used in.
        
        modeSymmetryX string {mustBeMember(modeSymmetryX, ...
            ["PEC", "PMC", "None"])} = "PMC";   % Symmetry of reflection about the x-axis.
        modeSymmetryY {mustBeTextScalar, mustBeMember(modeSymmetryY, ...
            ["PEC", "PMC", "None"])} = "PEC";   % Symmetry of reflection about the y-axis.
        modeSymmetryAxial {mustBeTextScalar, mustBeMember(modeSymmetryAxial, ...
            ["TE", "TM", "None"])} = "None";    % Axial symmetry.
    end
    properties (Dependent, Access=public)
        waveguideEr(1, :) {nLayer.mustBeErUrCallable};  % Array of filled waveguide permittivity values for each mode.
        waveguideUr(1, :) {nLayer.mustBeErUrCallable};  % Array of filled waveguide permeability values for each mode.

        excitationModeIndices(1, :) {mustBeInteger, mustBePositive};    % Array of indices of excitation modes (i.e., 'n' in Smn).
        receiveModeIndices(1, :) {mustBeInteger, mustBePositive};       % Array of indices of receiving modes (i.e., 'm' in Smn).
    end
    properties(Dependent, GetAccess=public, SetAccess=immutable)
        cutoffWavenumbers;  % Array of free-space cutoff wavenumbers (kc0) for each mode.
        modeTypes;          % String array of mode type for each mode ("TE", "TM", or "Hybrid").
        modeLabels;         % String array of labels for each mode.

        numModes;           % Total number of modes considered (TE + TM + Hybrid).
        numModes_TE;        % Number of TE modes considered.
        numModes_TM;        % Number of TM modes considered.
        numModes_Hybrid;    % Number of Hybrid modes considered.
    end
    properties (Access=private, Hidden)
        fixed_kr;       % Fixed-point integral coordindates kr.
        fixed_Ah;       % Fixed-point integral weights for Ah.
        fixed_Ae;       % Fixed-point integral weights for Ae.

        shouldRecomputeWeights(1, 1) logical = true;    % Flag to recompute integral weights.
    end

    %% Class Functions
    methods (Access=public)
        [outputLabels] = getOutputLabels(O);
    end
    methods (Access=protected)
        [modeStruct] = defineWaveguideModes(O);
    end
    methods (Access=protected)
        [gam] = calculate_impl(O, f, er, ur, thk);
    end
    methods (Access=private)
        [] = computeIntegralWeights(O, options);
        [krc, AhHat, AeHat] = computeAhat(O);
        [A] = computeA(O, f, er, ur, thk);
        [K] = computeK(O, f);
    end

    %% Class Constructor
    methods
        function O = nLayerOpenEnded(modeStructs)
            %NLAYEROPENENDED Construct an instance of this class.
            arguments (Repeating)
                modeStructs;
            end

            % % Pass in modeStructs
            % if isempty(modeStructs)
            %     error("Must pass in at least one modeStruct " + ...
            %         "to the constructor of 'nLayerOpenEnded'.");
            % end
            if ~isempty(modeStructs)
                O.modeStructs = modeStructs;
            end

        end
    end

    %% Class Setters
    methods
        function set.modeStructs(O, newModeStructs)
            O.modeStructs = newModeStructs;
            O.shouldRecomputeWeights = true;        %#ok<MCSUP>
            fprintf("Changed: modeStructs\n");
        end
        function set.frequencyRange(O, newFreqRange)
            O.frequencyRange = [min(newFreqRange); max(newFreqRange)];
            O.shouldRecomputeWeights = true;        %#ok<MCSUP>
            fprintf("Changed: freq\n");
        end

        function set.waveguideEr(O, newEr)
            if numel(newEr) ~= O.numModes && numel(newEr) ~= 1
                error("Array size must be equal to 1 or the number of modes.");
            end
            if ~isnumeric(newEr)
                newEr = {newEr};
            else
                newEr = mat2cell(newEr, ones(numel(newEr), 1));
            end
            if numel(newEr) == 1
                newEr = repmat(newEr, O.numModes, 1);
            end
            for ii = 1:numel(O.modeStructs)
                O.modeStructs(ii).WaveguideEr = newEr{ii};
            end
        end
        function set.waveguideUr(O, newUr)
            if numel(newUr) ~= O.numModes && numel(newUr) ~= 1
                error("Array size must be equal to 1 or the number of modes.");
            end
            if ~isnumeric(newUr)
                newUr = {newUr};
            else
                newUr = mat2cell(newUr, ones(numel(newUr), 1));
            end
            if numel(newUr) == 1
                newUr = repmat(newUr, O.numModes, 1);
            end
            for ii = 1:numel(O.modeStructs)
                O.modeStructs(ii).WaveguideUr = newUr{ii};
            end
        end
        function set.excitationModeIndices(O, newInds)
            isExMode = false(O.numModes, 1);
            isExMode(newInds) = true;
            [O.modeStructs.IsExcitationMode] = isExMode;
        end
        function set.receiveModeIndices(O, newInds)
            isRxMode = false(O.numModes, 1);
            isRxMode(newInds) = true;
            [O.modeStructs.IsReceiveMode] = isRxMode;
        end
    end

    %% Class Getters
    methods
        function [er] = get.waveguideEr(O)
            er = {O.modeStructs.WaveguideEr};
            if isequal(er{:})
                er = er{1};
            end
        end
        function [ur] = get.waveguideUr(O)
            ur = {O.modeStructs.WaveguideUr};
            if isequal(ur{:})
                ur = ur{1};
            end
        end
        function [inds] = get.excitationModeIndices(O)
            inds = find([O.modeStructs.IsExcitationMode]);
        end
        function [inds] = get.receiveModeIndices(O)
            inds = find([O.modeStructs.IsReceiveMode]);
        end

        function [kc0] = get.cutoffWavenumbers(O)
            kc0 = [O.modeStructs.CutoffWavenumber];
        end
        function [types] = get.modeTypes(O)
            types = [O.modeStructs.ModeType];
        end
        function [labels] = get.modeLabels(O)
            labels = [O.modeStructs.ModeLabel];
        end
        function [num] = get.numModes(O)
            num = numel(O.modeStructs);
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

