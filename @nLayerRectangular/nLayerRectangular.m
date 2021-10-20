classdef nLayerRectangular < nLayerForward
    %NLAYERRECTANGULAR Implementation of nLayerForward for rectangular waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by a rectangular waveguide looking into a multilayer structure.
    %   
    % Example Usages:
    %   NL = nLayerRectangular(band, maxM, maxN);
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, ur, thk, AbsTol);
    %
    %   NL = nLayerRectangular(a, b, maxM, maxN);
    %   gam = NL.calculate(f, er, ur, thk);
    %
    % There are a number of parameters that can be modified to change the
    % behavior. However, after changing any of the following parameters,
    % the member function "recomputeInterpolants()" must be called.
    %
    % List of critical parameters:
    %   a;
    %   b;
    %   c; (Inherited from nLayerForward)
    %   modesTE;
    %   interpolationPointsTau;
    %   integralPointsTauFixed;
    %   integralPointsPsi;
    %
    % It should be noted that changing the speed of light ("c") changes
    % the units used for all calculations, and thus "a" and "b" should be
    % changed as well.
    %
    % Author: Matt Dvorsky
    
    properties (GetAccess = public, SetAccess = public)
        a;          % Waveguide broad dimension (mm default units).
        b;          % Waveguide narrow dimension (mm default units).
        
        modesTE;    % List of TE modes to consider (in rows of m, n). First row must be [1, 0].
        
        interpolationPointsTau = 2^12;      % Number of points to use for the interpolation function along tau.
        integralPointsTauFixed = 135;       % Number of points to use for fixed-point integration along tau.
        integralInitialSegmentCount = 9;    % Initial number of segments along tau used in the adaptive integration routine.
        integralPointsPsi = 50;             % Number of points to use to integrate along psi.

        convergenceAbsTol = 0.001;  % Default tolerance for reflection coefficient calculations.
    end
    properties (GetAccess = public, SetAccess = private)
        numModes;   % Number of modes (numTE + numTM) to consider.
        modesTM;    % List of TM modes to consider (in rows of m, n).
    end
    properties (Access = private)
        integralScaleFactor;    % Scale factor for change of varibles from tau [0, inf) to tauP [0, 1].
        
        A1_EH;
        
        fixed_tau;
        fixed_A1_E;
        fixed_A1_H;
        fixed_errA1_E;
        fixed_errA1_H;
        
        init_tau;
        init_A1_E;
        init_A1_H;
                
        A2;
        b2;
    end
    
    %% Public member function definitions (implemented in separate files)
    methods (Access = public)
        gam = calculate(O, f, er, ur, thk, options);
        
        [modes] = setModes(O, maxM, maxN);
        [a, b] = setWaveguideBand(O, band);
        setWaveguideDimensions(O, a, b);
        recomputeInterpolants(O);
    end
    
    %% Private member function definitions (implemented in separate files)
    methods (Access = private)
        [A1_EH] = integrandA1(O, tauP, k0, er, ur, thk);
        [A1, b1] = computeA1b1(O, f, er, ur, thk, AbsTol)
        [tau, weightsE, weightsH] = computeIntegralInterpolant(O, orderTau);
        [A1, A2, b1, b2] = constructMatrixEquation(O, nLayerInt);
        [k_A1, k_A2, k_b1, k_b2] = constructFrequencyMultipliers(O, f);
    end
    
    %% Private static function definitions (implemented in separate files)
    methods (Static, Access = public)
        [specE, specH] = multilayerSpectrumRect(tau, k0, er, ur, thk);
    end
    
    %% Class constructor
    methods
        function O = nLayerRectangular(varargin)
            %NLAYERRECTANGULAR Construct an instance of this class
            % Usage:
            %   nLayerRectangular(a, b, maxM, maxN);
            %   nLayerRectangular(band, maxM, maxN);
            %
            % When calling as nLayerRectangular(a, b, maxM, maxN), the 
            % inputs a and b are the broad and narrow dimensions,
            % respectively, in mm, and maxM and maxN are the highest
            % mode indices to consider.
            %
            % When calling as nLayerRectangular(band, maxM, maxN), band
            % is a case-insensitive waveguide band identifier, and and 
            % maxM and maxN are the highest mode indices to consider.
            
            switch nargin
                case 3
                    O.setWaveguideBand(varargin{1});
                    O.setModes(varargin{2}, varargin{3});
                case 4
                    O.setWaveguideDimensions(varargin{1}, varargin{2});
                    O.setModes(varargin{3}, varargin{4});
                otherwise
                    error(strcat("The constructor of ", ...
                        "nLayerRectangular should be called using ", ...
                        "nLayerRectangular(a, b, maxM, maxN) ", ...
                        "or nLayerRectangular(band, maxM, maxN)."));
            end
            
            O.recomputeInterpolants();
        end
    end

end

