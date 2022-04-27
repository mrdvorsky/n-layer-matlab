classdef nLayerCircularTE < nLayerForward
    %NLAYERCIRCULARTE Implementation of nLayerForward for circular waveguides TE0n modes.
    % This class can be used to calculate the reflection coefficient seen
    % by a circular waveguide looking into a multilayer structure excited
    % with a TE01 mode.
    %
    % Example Usage:
    %   NL = nLayerCircularTE(numModes, R=wg_r);
    %   NL = nLayerCircularTE(numModes, R=wg_r, Verbosity=1);
    %   NL = nLayerCircularTE(numModes, R=wg_r, ...
    %           ConvergenceAbsTol=1e-4, IntegralPointsTauFixed=500);
    %   NL = nLayerCircularTE(numModes, R=wg_r, Prop=val, ...);
    %
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, [], thk);
    %   gam = NL.calculate(f, [], ur, thk);
    %   gam = NL.calculate(f, er, [], thk, BackingConductivity=sigma);
    %
    % nLayerRectangular Properties:
    %   r - Waveguide radius (mm default units).
    %   modesTE - List of TE modes to consider (in rows of m, n). First row
    %       must be [1, 0].
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
    %   integralInitialSegmentCount (9) - Initial number of segments along
    %       tau used in the adaptive integration routine. Must be an odd
    %       integer.
    %
    % There are a number of parameters that can be modified to change the
    % behavior. However, after changing any of the following parameters,
    % the member function "recomputeInterpolants()" must be called. This
    % function is automatically called upon construction of the object.
    % For Example:
    %   NL = nLayerCircularTE(numModes, R=wgR);
    %   NL.modesTE = [1, 0; 1, 2; 3, 0; 3, 2];
    %   NL.integralPointsTauFixed = 100;
    %   NL.recomputeInterpolants();     % Necessary in this case.
    %   [...] = NL.calculate(...);
    %
    % List of the critical parameters referenced above:
    %   r;
    %   speedOfLight; (Inherited from nLayerForward)
    %   modesTE;
    %   interpolationPointsTau;
    %   integralPointsTauFixed;
    %   integralInitialSegmentCount;
    %   integralPointsPsi;
    %
    % Any of the above properties can also be directly specified in the
    % class constructor: NL = nLayerRectangular(..., Prop=val, ...).
    % Constructing using this method avoids the requirement of having to
    % call recomputeInterpolants() manually.
    %
    % It should be noted that changing the speed of light ("speedOfLight")
    % changes the units used for all calculations, and thus "r" should be
    % changed as well.
    %
    % Author: Matt Dvorsky
    
    properties (GetAccess = public, SetAccess = public)
        r;                              % Waveguide radius (mm).
        modesTE;                        % List of TE0n modes (vector of n).
        interpolationPointsTau = 2^12;  % Number of points for lookup table along tau.
        integralPointsTauFixed = 300;   % Number of points for fixed point integral along tau.
        integralInitialSegmentCount = 9;    % Number of segments to start with in adaptive integral.
        integralPointsPsi = 50;         % Number of points for fixed point integral along psi.
        convergenceAbsTol = 0.001;      % Convergence tolerance value (absolute).
    end
    properties (GetAccess = public, SetAccess = private)
        numModes;           % Number of modes considered (of form TE0n).
        modeCutoffs;        % Cutoff wavenumbers of TE0n modes.
    end
    properties (Access = private)
        integralScaleFactor;    % Scale factor for change of varibles from 
                                % tau [0, inf) to tauP [0, 1].
        
        A1_H;           % Interpolation functions for A1_H(tauP).
        
        fixed_tau;      % Fixed-point integral coordindates tau.
        fixed_A1_H;     % Fixed-point integral weights for A1_H(tauP).
        fixed_errA1_H;  % Fixed-point error weights for A1_H(tauP).
        
        init_tau;       % First pass integral coordindates tau.
        init_A1_H;      % First pass preinterpolated A1_H(tauP).
        
        A2;             % Mode excitation matrix.
    end
    
    %% Protected member function definitions (implemented in separate files)
    methods (Access = protected)
        [gam] = calculateGamma(O, f, er, ur, thk);
    end
    
    %% Public member function definitions (implemented in separate files)
    methods (Access = public)
        [outputLabels] = getOutputLabels(O);
        
        setWaveguideDimensions(O, r);
        recomputeInterpolants(O);
    end
    
    %% Private member function definitions (implemented in separate files)
    methods (Access = private)
        [A1] = computeA1(O, f, er, ur, thk);
        [k_A1, k_A2, k_b1, k_b2] = constructFrequencyMultipliers(O, f);
        [A1_H] = integrandA1(O, tauP, k0, er, ur, thk);
        [integrandH] = computeIntegrandH(O, tauP);
        [A1, A2] = constructMatrixEquation(O, nLayerInt);
    end
    
    %% Private static function definitions (implemented in separate files)
    methods (Static, Access = public)
        [specH] = multilayerSpectrumH(tau, k0, er, ur, thk);
    end
    
    %% Class constructor
    methods
        function O = nLayerCircularTE(numModes, options)
            %NLAYERCIRCULARTE Construct an instance of this class.
            % Example Usage:
            %   NL = nLayerCircularTE(numModes, R=wg_r);
            %   NL = nLayerCircularTE(numModes, R=wg_r, Verbosity=1);
            %   NL = nLayerCircularTE(numModes, R=wg_r, ...
            %           ConvergenceAbsTol=1e-4, IntegralPointsTauFixed=500);
            %   NL = nLayerCircularTE(numModes, R=wg_r, Prop=val, ...);
            %
            % Inputs:
            %   numModes - Number of TE0n modes to consider.
            % Named Arguments:
            %   R (1) - Waveguide radius (mm).
            %   ModesTE - List of TE0n mode indices to use (vector of n).
            %       If specified, this will be used instead of numModes.
            %       The first value must be 1.
            %   Verbosity - Verbosity level. Set to 1 for basic command line
            %       output. Set to 2 for a per-frequency output.
            %   IntegralPointsTauFixed (300) - Number of points to use for
            %       fixed-point integration along tau.
            %   InterpolationPointsTau (2^12) - Number of points to use for the
            %       interpolation function along tau.
            %   IntegralInitialSegmentCount (9) - Initial number of
            %       segments along tau used in the adaptive integration
            %       routine. Must be an odd integer.
            %   IntegralPointsPsi (50) - Number of points to use to
            %       integrate along psi.
            %   ConvergenceAbsTol (0.001) - Default tolerance for
            %       reflection coefficient calculations.
            
            arguments
                numModes(1, 1) {mustBeInteger, mustBePositive};
                options.R(1, 1) {mustBeNumeric, mustBePositive} = 1;
                options.ModesTE(:, 2) {mustBeInteger, mustPositive};
                options.SpeedOfLight(1, 1) {mustBeReal, mustBePositive} = 299.792458;
                options.Verbosity(1, 1) {mustBeNumeric, mustBeNonnegative} = 0;
                options.ConvergenceAbsTol(1, 1) {mustBeNumeric, mustBePositive} = 0.001;
                options.IntegralPointsPsi(1, 1) {mustBeInteger, mustBePositive} = 50;
                options.IntegralPointsTauFixed(1, 1) {mustBeInteger, mustBePositive} = 300;
                options.InterpolationPointsTau(1, 1) {mustBeInteger, mustBePositive} = 2^12;
                options.IntegralInitialSegmentCount(1, 1) {mustBeInteger, mustBePositive} = 9;
            end
            
            %% Set Class Parameter Values
            if isfield(options, "modesTE")
                O.modesTE = options.modesTE;
            else
                O.modesTE = 1:numModes;
            end
            
            O.r = options.R;
            O.speedOfLight = options.SpeedOfLight;
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

