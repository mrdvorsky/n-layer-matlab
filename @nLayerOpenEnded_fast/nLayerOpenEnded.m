classdef nLayerOpenEnded < nLayerForward
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
    %   convergenceAbsTol (0.001) - Default tolerance for reflection
    %       coefficient calculations.
    %   integralPointsFixed_kRho (300) - Number of points to use for
    %       fixed-point integration along kRho. This parameter default is
    %       tuned so that the fixed-point method is used for loss tangents
    %       above ~0.1. Raise to lower this loss tangent threshold.
    %   interpolationPoints_kRho (2^12) - Number of points to use for the
    %       interpolation lookup table along kRho.
    %   integralPoints_kPhi (50) - Number of points to use to integrate
    %       along kPhi.
    %   integralInitialSegmentCount (9) - Initial number of segments along
    %       kRho used in the adaptive integration routine. Must be an odd
    %       integer.
    %   checkStructureValues (true) - Flag used in the "verifyStructure"
    %       function. If true, this function will throw errors if
    %       non-physical values of er, ur, or thk are passed in.
    %
    % Author: Matt Dvorsky

    properties (GetAccess=public, SetAccess=public)
        waveguideEr(1, 1) = 1;
        waveguideUr(1, 1) = 1;

        outputModes_TE(:, 1) logical = [];
        outputModes_TM(:, 1) logical = [];
        outputModes_Hybrid(:, 1) logical = [];

        interpolationPoints_kRho(1, 1) {mustBePositive, mustBeInteger} = 2^12;  % Number of points for lookup table along kRho.
        integralPointsFixed_kRho(1, 1) {mustBePositive, mustBeInteger} = 300;   % Number of points for fixed point integral along kRho.
        integralPoints_kPhi(1, 1) {mustBePositive, mustBeInteger} = 50;         % Number of points for fixed point integral along kPhi.
        integralInitialSegmentCount(1, 1) {mustBePositive, mustBeInteger} = 9;  % Number of segments to start with in adaptive integral.
        convergenceAbsTol(1, 1) {mustBePositive} = 0.001;                       % Convergence tolerance value.
    end
    properties (GetAccess=public, SetAccess=protected)
        numModes;           % Total number of modes considered (TE + TM + Hybrid).
        numModes_TE;        % Number of TE modes considered.
        numModes_TM;        % Number of TM modes considered.
        numModes_Hybrid;    % Number of Hybrid modes considered.

        cutoffBeta_TE(:, 1) = [];
        cutoffBeta_TM(:, 1) = [];
        cutoffBeta_Hybrid(:, 1) = [];
    end
    properties (Access=private, Hidden)
        integralScaleFactor = 1;    % Scale factor for change of varibles from kRho [0, inf) to kRhoP [0, 1].

        table_AheHat;       % Interpolation tables for AhHat(kRhoP) and AeHat(kRhoP).

        fixed_kRho;         % Fixed-point integral coordindates kRho.
        fixed_AhHat;        % Fixed-point integral weights for AhHat(kRhoP).
        fixed_AeHat;        % Fixed-point integral weights for AeHat(kRhoP).
        fixed_errorAhHat;   % Fixed-point error weights for AhHat(kRhoP).
        fixed_errorAeHat;   % Fixed-point error weights for AeHat(kRhoP).

        init_kRho;      % First pass integral coordindates kRho.
        init_AhHat;     % First pass preinterpolated AhHat(kRhoP).
        init_AeHat;     % First pass preinterpolated AeHat(kRhoP).
    end

    %% Class Functions
    methods (Access=public)
        [outputLabels] = getOutputLabels(O);
        recomputeInterpolants(O, options);
    end
    methods (Static, Access=public)
        [Ahat] = integrandAhat(O, kRhoP, k0, er, ur, thk);
        [Gamma0h, Gamma0e] = computeGamma0(kRho, k0, er, ur, thk);

        [modeStruct] = createModeStruct(options);
        [] = plotModeStruct(modeStruct, options);

        [modeStruct] = createModeStructRectangular(options);
        [modeStruct] = createModeStructCircular(options);
        [modeStruct] = createModeStructCoaxial(options);

        [Ex, Ey, cutoffWavenumber, phaseScale] = getSpectrumRectangular(wgA, wgB, m, n, TE_TM);
        [Ex, Ey, cutoffWavenumber, phaseScale] = getSpectrumCircular(wgR, m, n, TE_TM);
        [Ex, Ey, cutoffWavenumber, phaseScale] = getSpectrumCoaxial(wgR_inner, wgR_outer, m, n, TE_TM);
    end

    methods (Access=protected)
        [modeStruct] = defineWaveguideModes(O);
    end

    methods (Access=protected)
        [gam] = calculate_impl(O, f, er, ur, thk);
    end

    methods (Access=private)
        [AhHat, AeHat] = computeAhat(O, kRhoP);
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

            %% Pass in modeStructs
            % if isempty(modeStructs)
            %     error("Must pass in at least one modeStruct " + ...
            %         "to the constructor of 'nLayerOpenEnded'.");
            % end
            O.recomputeInterpolants(modeStructs{:});

        end
    end

end

