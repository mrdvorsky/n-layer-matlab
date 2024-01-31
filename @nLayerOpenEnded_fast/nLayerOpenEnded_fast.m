classdef nLayerOpenEnded_fast < nLayerForward
    %NLAYEROPENENDED Implementation of nLayerForward for general open-ended waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by an arbitrary waveguide looking into a multilayer structure. This
    % class is meant to be subclassed, where the subclass defines the modes
    % to be considered.
    %
    % Note that the units of all parameters should match that of the speed
    % of light specified by the speedOfLight parameter (default is mm GHz).
    %
    % nLayerOpenEnded Properties:
    %   speedOfLight (299.792458) - Speed of light (default is mm GHz).
    %   verbosity (0) - Verbosity level. Set to 1 for basic command line
    %       output. Set to 2 for a per-frequency output.
    %
    % Author: Matt Dvorsky

    properties (GetAccess=public, SetAccess=public)
        waveguideEr(1, 1) = 1;
        waveguideUr(1, 1) = 1;

        outputModes_TE(:, 1) logical = [];
        outputModes_TM(:, 1) logical = [];
        outputModes_Hybrid(:, 1) logical = [];
    end
    properties (GetAccess=public, SetAccess=protected)
        numModes;           % Total number of modes considered (TE + TM + Hybrid).
        numModes_TE;        % Number of TE modes considered.
        numModes_TM;        % Number of TM modes considered.
        numModes_Hybrid;    % Number of Hybrid modes considered.

        modeStructs;

        cutoffBeta_TE(:, 1) = [];
        cutoffBeta_TM(:, 1) = [];
        cutoffBeta_Hybrid(:, 1) = [];
    end
    properties (Access=private, Hidden)
        fixed_kRho;             % Fixed-point integral coordindates kRho.
        fixed_Ah_weights;       % Fixed-point integral weights for Ah.
        fixed_Ae_weights;       % Fixed-point integral weights for Ae.
    end

    %% Class Functions
    methods (Access=public)
        [outputLabels] = getOutputLabels(O);
        recomputeInterpolants(O, options);
    end

    methods (Access=protected)
        [modeStruct] = defineWaveguideModes(O);
    end

    methods (Access=protected)
        [gam] = calculate_impl(O, f, er, ur, thk);
    end

    methods (Access=private)
        [krc, AhHat, AeHat] = computeAhat(O, modeStruct);
        [A] = computeA(O, f, er, ur, thk);
        [K] = computeK(O, f);
    end

    %% Class Constructor
    methods
        function O = nLayerOpenEnded_fast(modeStructs)
            %NLAYEROPENENDED Construct an instance of this class.

            arguments (Repeating)
                modeStructs;
            end

            %% Pass in modeStructs
            % if isempty(modeStructs)
            %     error("Must pass in at least one modeStruct " + ...
            %         "to the constructor of 'nLayerOpenEnded'.");
            % end
            if ~isempty(modeStructs)
                O.modeStructs = modeStructs;
                O.recomputeInterpolants();
            end

        end
    end

end

