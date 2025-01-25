classdef nLayerCircularTE_old < nLayerForward
    %NLAYERCIRCULARTE Implementation of nLayerForward for circular waveguides TE0n modes.
    % This class can be used to calculate the reflection coefficient seen
    % by a circular waveguide looking into a multilayer structure, in
    % addition to the full mode S-parameter matrix. Note that the units of
    % all parameters should match that of the speed of light specified by
    % the speedOfLight parameter (default is mm GHz).
    %
    % Example Usage:
    %   NL = nLayerCircularTE(numModes, waveguideR=5.8);
    %   NL = nLayerCircularTE(numModes, waveguideR=5.8, verbosity=1);
    %   NL = nLayerCircularTE(numModes, waveguideR=5.8e-3, ...
    %       speedOfLight=299.79e6);
    %   NL = nLayerCircularTE(maxM, maxN, waveguideBand="ka", ...
    %       prop1=val1, prop2=val2, ...);
    %
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, [], thk);
    %   gam = NL.calculate(f, [], ur, thk);
    %   gam = NL.calculate(f, er, [], thk, BackingConductivity=sigma);
    %
    % nLayerCircularTE Properties:
    %   waveguideR - Waveguide radius.
    %   speedOfLight (299.792458) - Speed of light (default is mm GHz).
    %   modesTE - List of TE0n modes to consider (in rows of [n]). The
    %       ports for the S-parameters will be in the order specified by
    %       modesTE.
    %   numModes (read-only) - Number of TE0n modes to consider.
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
    %   NL = nLayerCircularTE(numModes, waveguideR=wgR);
    %   NL.modesTE = [1; 2; 3; 4];
    %   NL.integralPointsFixed_kRho = 100;
    %   NL.recomputeInterpolants(); % This line is necessary in this case.
    %   [...] = NL.calculate(...);
    %
    % List of the critical parameters referenced above:
    %   waveguideR;
    %   speedOfLight;
    %   modesTE;
    %   interpolationPoints_kRho;
    %   integralPointsFixed_kRho;
    %   integralInitialSegmentCount;
    %
    % Any of the above properties can also be directly specified in the
    % class constructor: NL = nLayerCircularTE(..., prop=val, ...).
    % Constructing using this method avoids the requirement of having to
    % call "recomputeInterpolants()" manually.
    %
    % Author: Matt Dvorsky

    properties (Access=public)
        waveguideR(1, 1) {mustBePositive} = 1;      % Waveguide radius.
        waveguideEr(1, 1) {nLayer.mustBeValidErUr} = 1;  % Waveguide fill er.
        waveguideUr(1, 1) {nLayer.mustBeValidErUr} = 1;  % Waveguide fill ur.
        modesTE(:, 1) {mustBeInteger, mustBeNonnegative};   % List of TE0n modes (vector of n).
        convergenceAbsTol(1, 1) {mustBePositive} = 0.001;   % Convergence tolerance value.

        interpolationPoints_kRho(1, 1) {mustBePositive, mustBeInteger} = 2^12;  % Number of points for lookup table along kRho.
        integralPointsFixed_kRho(1, 1) {mustBePositive, mustBeInteger} = 301;   % Number of points for fixed point integral along kRho.
        integralInitialSegmentCount(1, 1) {nLayer.mustBePositiveOddInteger} = 9; % Number of segments to start with in adaptive integral.
    end
    properties (GetAccess=public, SetAccess=private)
        numModes;           % Number of modes considered (of form TE0n).
        modeCutoffs;        % Cutoff wavenumbers of TE0n modes.
    end
    properties (Access=private)
        integralScaleFactor;    % Scale factor for change of varibles from kRho [0, inf) to kRhoP [0, 1].

        table_AhHat;        % Interpolation table for AhHat(kRhoP).

        fixed_kRho;         % Fixed-point integral coordindates kRho.
        fixed_AhHat;        % Fixed-point integral weights for AhHat(kRhoP).
        fixed_errorAhHat;   % Fixed-point error weights for AhHat(kRhoP).

        init_kRho;      % First pass integral coordindates kRho.
        init_AhHat;     % First pass preinterpolated AhHat(kRhoP).
    end

    %% Class Functions
    methods (Access=public)
        [gam] = calculate(O, f, er, ur, thk);
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

        [AhHat] = computeAhat(O, kRhoP);
        [Ahat] = integrandAhat(O, kRhoP, k0, er, ur, thk);
    end
    methods (Static, Access=public)
        [Gamma0h] = computeGamma0(kRho, k0, er, ur, thk);
    end

    %% Class constructor
    methods
        function O = nLayerCircularTE_old(numModes, classProperties)
            %NLAYERCIRCULARTE Construct an instance of this class.
            % Example Usage:
            %   See example usage in main class documentation. Note that
            %   all public class properties can be specified as a named
            %   argument to the constructor (e.g., as "verbosity=1").
            %
            % Inputs:
            %   numModes - Number of TE0n modes to consider.

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
                O.(classProperties{ii}) = classProperties{ii + 1};
            end

            if isempty(O.modesTE)
                O.modesTE = (1:numModes).';
            end

            O.recomputeInterpolants();
        end
    end

end


