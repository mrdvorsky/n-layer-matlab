classdef nLayerCircularTM < nLayerForward
    %NLAYERCIRCULARTM Implementation of nLayerForward for circular waveguides TM0n modes.
    % This class can be used to calculate the reflection coefficient(s)
    % seen by a circular waveguide looking into a multilayer structure, in
    % addition to the full mode S-parameter matrix. Note that the units of
    % all parameters should match that of the speed of light specified by
    % the speedOfLight parameter (default is mm GHz).
    %
    % Example Usage:
    %   NL = nLayerCircularTM(numModes, waveguideR=4.5);
    %   NL = nLayerCircularTM(numModes, waveguideR=4.5, verbosity=1);
    %   NL = nLayerCircularTM(numModes, waveguideR=4.5e-3, ...
    %       speedOfLight=299.79e6);
    %   NL = nLayerCircularTM(maxM, maxN, waveguideBand="ka", ...
    %       prop1=val1, prop2=val2, ...);
    %
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, [], thk);
    %   gam = NL.calculate(f, [], ur, thk);
    %   gam = NL.calculate(f, er, [], thk, BackingConductivity=sigma);
    %
    % nLayerCircularTM Properties:
    %   waveguideR - Waveguide radius.
    %   speedOfLight (299.792458) - Speed of light (default is mm GHz).
    %   modesTM - List of TM0n modes to consider (in rows of [n]). The
    %       ports for the S-parameters will be in the order specified by
    %       modesTM.
    %   numModes (read-only) - Number of TM0n modes to consider.
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
    %   integralInitialSegmentCount (9) - Initial number of segments along
    %       kRho used in the adaptive integration routine. Must be an odd
    %       integer.
    %   checkStructureValues (true) - Flag used in the "verifyStructure"
    %       function. If true, this function will throw errors if
    %       non-physical values of er, ur, or thk are passed in.
    %
    % There are a number of parameters that can be modified to change the
    % behavior. However, after changing any of the following parameters,
    % the member function "recomputeInterpolants()" must be called. This
    % function is automatically called upon construction of the object.
    % For Example:
    %   NL = nLayerCircularTM(numModes, waveguideR=wgR);
    %   NL.modesTM = [1; 2; 3; 4];
    %   NL.integralPointsFixed_kRho = 100;
    %   NL.recomputeInterpolants(); % This line is necessary in this case.
    %   [...] = NL.calculate(...);
    %
    % List of the critical parameters referenced above:
    %   waveguideR;
    %   speedOfLight;
    %   modesTM;
    %   interpolationPoints_kRho;
    %   integralPointsFixed_kRho;
    %   integralInitialSegmentCount;
    %
    % Any of the above properties can also be directly specified in the
    % class constructor: NL = nLayerCircularTM(..., prop=val, ...).
    % Constructing using this method avoids the requirement of having to
    % call "recomputeInterpolants()" manually.
    %
    % Author: Matt Dvorsky

    properties (Access=public)
        waveguideR(1, 1) {mustBePositive} = 1;      % Waveguide radius.
        waveguideEr(1, 1) {nLayerForward.mustBeValidErUr} = 1;  % Waveguide fill er.
        waveguideUr(1, 1) {nLayerForward.mustBeValidErUr} = 1;  % Waveguide fill ur.
        modesTM(:, 1) {mustBeInteger, mustBeNonnegative};   % List of TM0n modes (vector of n).
        convergenceAbsTol(1, 1) {mustBePositive} = 0.001;   % Convergence tolerance value.

        interpolationPoints_kRho(1, 1) {mustBePositive, mustBeInteger} = 2^12;  % Number of points for lookup table along kRho.
        integralPointsFixed_kRho(1, 1) {mustBePositive, mustBeInteger} = 300;   % Number of points for fixed point integral along kRho.
        integralInitialSegmentCount(1, 1) {nLayerForward.mustBePositiveOddInteger} = 9; % Number of segments to start with in adaptive integral.
    end
    properties (GetAccess=public, SetAccess=private)
        numModes;           % Number of modes considered (of form TM0n).
        modeCutoffs;        % Cutoff wavenumbers of TM0n modes.
    end
    properties (Access=private)
        integralScaleFactor;    % Scale factor for change of varibles from kRho [0, inf) to kRhoP [0, 1].

        table_AeHat;        % Interpolation table for AeHat(kRhoP).

        fixed_kRho;         % Fixed-point integral coordindates kRho.
        fixed_AeHat;        % Fixed-point integral weights for AeHat(kRhoP).
        fixed_errorAeHat;   % Fixed-point error weights for AeHat(kRhoP).

        init_kRho;      % First pass integral coordindates kRho.
        init_AeHat;     % First pass preinterpolated AeHat(kRhoP).
    end

    %% Class Functions
    methods (Access=protected)
        [gam] = calculate_impl(O, f, er, ur, thk);
    end
    methods (Access=public)
        [outputLabels] = getOutputLabels(O);
        setWaveguideDimensions(O, waveguideA, waveguideB);
        recomputeInterpolants(O);
    end
    methods (Access=private)
        [A] = computeA(O, f, er, ur, thk);
        [B] = computeB(O);
        [kA, kB] = computeK(O, f);

        [AeHat] = computeAhat(O, kRhoP);
        [Ahat] = integrandAhat(O, kRhoP, k0, er, ur, thk);
    end
    methods (Static, Access=public)
        [Gamma0e] = computeGamma0(kRho, k0, er, ur, thk);
    end

    %% Class constructor
    methods
        function O = nLayerCircularTM(numModes, classProperties)
            %NLAYERCIRCULARTM Construct an instance of this class.
            % Example Usage:
            %   See example usage in main class documentation. Note that
            %   all public class properties can be specified as a named
            %   argument to the constructor (e.g., as "verbosity=1").
            %
            % Inputs:
            %   numModes - Number of TM0n modes to consider.

            arguments
                numModes(1, 1) {mustBeInteger, mustBePositive};
            end
            arguments (Repeating)
                classProperties;
            end

            %% Set Class Parameter Values
            if mod(numel(classProperties), 2) ~= 0
                error("Parameter and value arguments must come in pairs.");
            end
            for ii = 1:2:numel(classProperties)
                set(O, classProperties{ii}, classProperties{ii + 1});
            end

            if isempty(O.modesTM)
                O.modesTM = (1:numModes).';
            end

            O.recomputeInterpolants();
        end
    end

end

