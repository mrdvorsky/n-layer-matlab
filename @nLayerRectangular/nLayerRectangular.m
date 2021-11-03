classdef nLayerRectangular < nLayerForward
    %NLAYERRECTANGULAR Implementation of nLayerForward for rectangular waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by a rectangular waveguide looking into a multilayer structure.
    %
    % Example Usage:
    %   NL = nLayerRectangular(maxM, maxN, A=7.112, B=3.556);
    %   NL = nLayerRectangular(maxM, maxN, Band="x");
    %   NL = nLayerRectangular(maxM, maxN, Band="ka", Verbosity=1);
    %   NL = nLayerRectangular(maxM, maxN, Band="v", ...
    %           ConvergenceAbsTol=1e-4, IntegralPointsTauFixed=500);
    %   NL = nLayerRectangular(maxM, maxN, band="x", prop=val, ...);
    %
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, ur, thk, AbsTol);
    %
    % nLayerRectangular Properties:
    %   a - Waveguide broad dimension (mm default units).
    %   b - Waveguide narrow dimension (mm default units).
    %   modesTE - List of TE modes to consider (in rows of m, n). First row
    %       must be [1, 0].
    %   modesTM (read-only) - List of TM modes to consider (in rows of m, n).
    %   numModes (read-only) - Number of modes (numTE + numTM) to consider.
    %   verbosity - Verbosity level. Set to 1 for basic command line
    %       output. Set to 2 for a per-frequency output.
    %   convergenceAbsTol (0.001) - Default tolerance for reflection
    %       coefficient calculations.
    %   integralPointsTauFixed (300) - Number of points to use for
    %       fixed-point integration along tau. This parameter default is
    %       tuned so that the fixed-point method is used for loss tangents
    %       above ~0.1. Raise to lower this loss tangent threshold.
    %   interpolationPointsTau (2^12) - Number of points to use for the
    %       interpolation function along tau.
    %   integralPointsPsi (50) - Number of points to use to integrate along psi.
    %   integralInitialSegmentCount (9) - Initial number of segments along
    %       tau used in the adaptive integration routine. Must be an odd
    %       integer.
    %
    % There are a number of parameters that can be modified to change the
    % behavior. However, after changing any of the following parameters,
    % the member function "recomputeInterpolants()" must be called. This
    % function is automatically called upon construction of the object.
    % For Example:
    %   NL = nLayerRectangular(maxM, maxN, Band="x");
    %   NL.modesTE = [1, 0; 1, 2; 3, 0; 3, 2];
    %   NL.integralPointsTauFixed = 100;
    %   NL.recomputeInterpolants();     % This line is necessary.
    %   [...] = NL.calculate(...);
    %
    % List of the critical parameters referenced above:
    %   a;
    %   b;
    %   c; (Inherited from nLayerForward)
    %   modesTE;
    %   interpolationPointsTau;
    %   integralPointsTauFixed;
    %   integralPointsPsi;
    %   integralInitialSegmentCount;
    %
    % Any of the above properties can also be directly specified in the
    % class constructor: NL = nLayerRectangular(..., Prop=val, ...).
    % Constructing using this method avoids the requirement of having to
    % call recomputeInterpolants() manually.
    %
    % It should be noted that changing the speed of light ("c") changes
    % the units used for all calculations, and thus "a" and "b" should be
    % changed as well.
    %
    % Author: Matt Dvorsky
    
    properties (GetAccess = public, SetAccess = public)
        a;
        b;
        modesTE;
        interpolationPointsTau = 2^12;
        integralPointsTauFixed = 300;
        integralInitialSegmentCount = 9;
        integralPointsPsi = 50;
        convergenceAbsTol = 0.001;
    end
    properties (GetAccess = public, SetAccess = private)
        numModes;
        modesTM;
    end
    properties (Access = private)
        integralScaleFactor;    % Scale factor for change of varibles from 
                                % tau [0, inf) to tauP [0, 1].
        
        A1_EH;  % Interpolation functions for A1_E(tauP) and A1_H(tauP).
        
        fixed_tau;      % Fixed-point integral coordindates tau.
        fixed_A1_E;     % Fixed-point integral weights for A1_E(tauP).
        fixed_A1_H;     % Fixed-point integral weights for A1_H(tauP).
        fixed_errA1_E;  % Fixed-point error weights for A1_E(tauP).
        fixed_errA1_H;  % Fixed-point error weights for A1_H(tauP).
        
        init_tau;       % First pass integral coordindates tau.
        init_A1_E;      % First pass preinterpolated A1_E(tauP).
        init_A1_H;      % First pass preinterpolated A1_H(tauP).
                
        A2;             % Mode excitation matrix.
        b2;             % Mode excitation vector. Equal to A2(:, 1, ...).
    end
    
    %% Public member function definitions (implemented in separate files)
    methods (Access = public)
        gam = calculate(O, f, er, ur, thk, options);
        
        [modes] = setModes(O, maxM, maxN);
        [a, b] = setWaveguideBand(O, band, options);
        setWaveguideDimensions(O, a, b);
        recomputeInterpolants(O);
    end
    
    %% Private member function definitions (implemented in separate files)
    methods (Access = private)
        [A1, b1] = computeA1b1(O, f, er, ur, thk, AbsTol);
        [k_A1, k_A2, k_b1, k_b2] = constructFrequencyMultipliers(O, f);
        [A1_EH] = integrandA1(O, tauP, k0, er, ur, thk);
        [integrandE, integrandH] = computeIntegrandEH(O, tauP);
        [A1, A2, b1, b2] = constructMatrixEquation(O, nLayerInt);
    end
    
    %% Private static function definitions (implemented in separate files)
    methods (Static, Access = public)
        [specE, specH] = multilayerSpectrumEH(tau, k0, er, ur, thk);
    end
    
    %% Class constructor
    methods
        function O = nLayerRectangular(maxM, maxN, options)
            %NLAYERRECTANGULAR Construct an instance of this class.
            % Example Usage:
            %   NL = nLayerRectangular(maxM, maxN, A=7.112, B=3.556);
            %   NL = nLayerRectangular(maxM, maxN, Band="x");
            %   NL = nLayerRectangular(maxM, maxN, Band="x", Verbosity=1);
            %   NL = nLayerRectangular(maxM, maxN, Band="x", ...
            %       ConvergenceAbsTol=1e-4, IntegralPointsTauFixed=500);
            %   NL = nLayerRectangular(maxM, maxN, Band="x", Prop=val, ...);
            %
            % Inputs:
            %   maxM - Highest index m of TEmn and TMmn modes to consider.
            %   maxN - Highest index n of TEmn and TMmn modes to consider.
            % Named Options:
            %   A (1) - Waveguide broad dimension (mm).
            %   B (0.5) - Waveguide narrow dimension (mm).
            %   Band - Case-insensitive waveguide band to use. Either
            %       specify this or the dimensions (a and b) directly.
            %   ModesTE - List of modes to use in rows of [m, n].
            %       If specified, this will be used instead of maxM and
            %       maxN. The first row must be [1, 0].
            %   Verbosity - Verbosity level. Set to 1 for basic command line
            %       output. Set to 2 for a per-frequency output.
            %   ConvergenceAbsTol (0.001) - Default tolerance for
            %       reflection coefficient calculations.
            %   IntegralPointsTauFixed (300) - Number of points to use for
            %       fixed-point integration along tau.
            %   InterpolationPointsTau (2^12) - Number of points to use for the
            %       interpolation function along tau.
            %   IntegralPointsPsi (50) - Number of points to use to
            %       integrate along psi.
            %   IntegralInitialSegmentCount (9) - Initial number of
            %       segments along tau used in the adaptive integration
            %       routine. Must be an odd integer.
            
            arguments
                maxM(1, 1) {mustBeInteger, mustBePositive};
                maxN(1, 1) {mustBeInteger, mustBeNonnegative};
                options.A(1, 1) {mustBeNumeric, mustBePositive} = 1;
                options.B(1, 1) {mustBeNumeric, mustBePositive} = 0.5;
                options.Band {mustBeTextScalar};
                options.ModesTE(:, 2) {mustBeInteger, mustBeNonnegative};
                options.Verbosity(1, 1) {mustBeNumeric, mustBeNonnegative} = 0;
                options.ConvergenceAbsTol(1, 1) {mustBeNumeric, mustBePositive} = 0.001;
                options.IntegralPointsPsi(1, 1) {mustBeInteger, mustBePositive} = 50;
                options.IntegralPointsTauFixed(1, 1) {mustBeInteger, mustBePositive} = 300;
                options.InterpolationPointsTau(1, 1) {mustBeInteger, mustBePositive} = 2^12;
                options.IntegralInitialSegmentCount(1, 1) {mustBeInteger, mustBePositive} = 9;
            end
            
            %% Set Class Parameter Values
            if isfield(options, "ModesTE")
                O.modesTE = options.ModesTE;
            else
                O.setModes(maxM, maxN);
            end
            
            if isfield(options, "Band")
                O.setWaveguideBand(options.Band);
            else
                O.setWaveguideDimensions(options.A, options.B);
            end
            
            O.verbosity = options.Verbosity;
            O.convergenceAbsTol = options.ConvergenceAbsTol;
            O.integralPointsPsi = options.IntegralPointsPsi;
            O.integralPointsTauFixed = options.IntegralPointsTauFixed;
            O.interpolationPointsTau = options.InterpolationPointsTau;
            O.integralInitialSegmentCount = options.IntegralInitialSegmentCount;
            
            O.recomputeInterpolants();
        end
    end

end

