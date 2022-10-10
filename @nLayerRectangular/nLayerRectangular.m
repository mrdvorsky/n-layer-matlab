classdef nLayerRectangular < nLayerForward
    %NLAYERRECTANGULAR Implementation of nLayerForward for rectangular waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by a rectangular waveguide looking into a multilayer structure. Note
    % that the units of all parameters should match that of the speed of
    % light specified by the speedOfLight parameter (default is mm GHz).
    %
    % Example Usage:
    %   NL = nLayerRectangular(maxM, maxN, waveguideBand="x");
    %   NL = nLayerRectangular(maxM, maxN, waveguideA=7.112, waveguideB=3.556);
    %   NL = nLayerRectangular(maxM, maxN, speedOfLight=299.79e6, ...
    %       waveguideA=7.112e-3, waveguideB=3.556e-3);
    %   NL = nLayerRectangular(maxM, maxN, waveguideBand="ka", ...
    %       prop1=val1, prop2=val2, ...);
    %
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, [], thk);
    %   gam = NL.calculate(f, [], ur, thk);
    %   gam = NL.calculate(f, er, [], thk, BackingConductivity=sigma);
    %
    % nLayerRectangular Properties:
    %   waveguideA - Waveguide broad dimension.
    %   waveguideB - Waveguide narrow dimension.
    %   speedOfLight (299.792458) - Speed of light (default is mm GHz).
    %   modesTE - List of TEmn modes to consider (in rows of [m, n]). The
    %       ports for the S-parameters will be in the order specified by
    %       modesTE, followed by modesTM.
    %   modesTM (read-only) - List of TMmn modes to consider (in rows of 
    %       [m, n]). This is automatically generated from modesTE. It will
    %       be in the same order after removing invalid TM modes.
    %   numModes (read-only) - Number of modes (numTE + numTM) to consider.
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
    % There are a number of parameters that can be modified to change the
    % behavior. However, after changing any of the following parameters,
    % the member function "recomputeInterpolants()" must be called. This
    % function is automatically called upon construction of the object.
    % For Example:
    %   NL = nLayerRectangular(maxM, maxN, Band="x");
    %   NL.modesTE = [1, 0; 1, 2; 3, 0; 3, 2];
    %   NL.integralPointsFixed_kRho = 100;
    %   NL.recomputeInterpolants(); % This line is necessary in this case.
    %   [...] = NL.calculate(...);
    %
    % List of the critical parameters referenced above:
    %   waveguideA;
    %   waveguideB;
    %   speedOfLight;
    %   modesTE;
    %   interpolationPoints_kRho;
    %   integralPointsFixed_kRho;
    %   integralPoints_kPhi;
    %   integralInitialSegmentCount;
    %
    % Any of the above properties can also be directly specified in the
    % class constructor: NL = nLayerRectangular(..., prop=val, ...).
    % Constructing using this method avoids the requirement of having to
    % call "recomputeInterpolants()" manually.
    %
    % Author: Matt Dvorsky

    properties (GetAccess=public, SetAccess=public)
        waveguideA(1, 1) {mustBePositive} = 1;              % Waveguide broad dimension.
        waveguideB(1, 1) {mustBePositive} = 0.5;            % Waveguide narrow dimension.
        modesTE(:, 2) {mustBeInteger, mustBeNonnegative};   % List of TEmn modes in rows of [m, n].
        interpolationPoints_kRho(1, 1) {mustBePositive, mustBeInteger} = 2^12;  % Number of points for lookup table along kRho.
        integralPointsFixed_kRho(1, 1) {mustBePositive, mustBeInteger} = 300;   % Number of points for fixed point integral along kRho.
        integralPoints_kPhi(1, 1) {mustBePositive, mustBeInteger} = 50;         % Number of points for fixed point integral along kPhi.
        integralInitialSegmentCount(1, 1) {mustBePositive, mustBeInteger} = 9;  % Number of segments to start with in adaptive integral.
        convergenceAbsTol(1, 1) {mustBePositive} = 0.001;                       % Convergence tolerance value.
    end
    properties (GetAccess=public, SetAccess=private)
        numModes;       % Number of modes considered (TE + TM).
        modesTM;        % List of TMmn modes in rows of [m, n].
        waveguideBand = "";     % Waveguide band identifier.
    end
    properties (Access=private)
        integralScaleFactor;    % Scale factor for change of varibles from kRho [0, inf) to kRhoP [0, 1].

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
    methods (Access=protected)
        [gam] = calculate_impl(O, f, er, ur, thk);
    end
    methods (Access=public)
        [outputLabels] = getOutputLabels(O);

        [modesTE, modesTM] = setModes(O, maxM, maxN);
        [waveguideA, waveguideB] = setWaveguideBand(O, band, options);

        recomputeInterpolants(O);
    end
    methods (Access=private)
        [A] = computeA(O, f, er, ur, thk);
        [B] = computeB(O);
        [kA, kB] = computeK(O, f);

        [AhHat, AeHat] = computeAhat(O, kRhoP);
        [Ahat] = integrandAhat(O, kRhoP, k0, er, ur, thk);
    end
    methods (Static, Access=public)
        [Gamma0h, Gamma0e] = computeGamma0(kRho, k0, er, ur, thk);
    end

    %% Class constructor
    methods
        function O = nLayerRectangular(maxM, maxN, classProperties)
            %NLAYERRECTANGULAR Construct an instance of this class.
            % Example Usage:
            %   See example usage in main class documentation. Note that
            %   all public class properties can be specified as a named
            %   argument to the constructor (e.g., as "verbosity=1").
            %
            % Inputs:
            %   maxM - Highest index m of TEmn and TMmn modes to consider.
            %   maxN - Highest index n of TEmn and TMmn modes to consider.

            arguments
                maxM(1, 1) {mustBeInteger, mustBePositive};
                maxN(1, 1) {mustBeInteger, mustBeNonnegative};
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

            if isempty(O.modesTE)
                O.setModes(maxM, maxN);
            end

            if strlength(O.waveguideBand) > 0
                O.setWaveguideBand(O.waveguideBand);
            end

            O.recomputeInterpolants();
        end
    end

end

