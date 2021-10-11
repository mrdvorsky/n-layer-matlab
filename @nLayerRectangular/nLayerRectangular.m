classdef nLayerRectangular < nLayerForward
    %NLAYERRECTANGULAR Implementation of nLayerForward for rectangular waveguides.
    %   This class can be used to calculate the reflection coefficient seen
    %   by a rectangular waveguide looking into a multilayer structure.
    %   
    %   Example Usages:
    %       NL = nLayerRectangular(band, maxM, maxN);
    %       gam = NL.calculate(f, er, ur, thk);
    %       gam = NL.calculate(f, er, ur, thk, AbsTol);
    %
    %       NL = nLayerRectangular(a, b, maxM, maxN);
    %       gam = NL.calculate(f, er, ur, thk);
    %
    %   There are a number of parameters that can be modified to change the
    %   behavior. However, after changing any of the following parameters,
    %   the member function "recomputeInterpolants()" must be called.
    %   List of these parameters: "a", "b", "c", "modes",
    %       "interpolationPointsTau", "integralPointsPsi".
    %   It should be noted that changing the speed of light ("c") changes
    %   the units used for all calculations, and thus "a" and "b" should be
    %   changed as well.
    
    properties (Access = private)
        integralScaleFactor;
        
        A1b1_EH;
        
        lossy_tau;
        lossy_A1b1_E;
        lossy_A1b1_H;
        lossy_errA1b1_E;
        lossy_errA1b1_H;
        
        init_A1b1_E;
        init_A1b1_H;
                
        A2;
        b2;
    end
    properties (GetAccess = public, SetAccess = public)
        a;          % Waveguide broad dimension (mm default units).
        b;          % Waveguide narrow dimension (mm default units).
        
        modes;      % List of TE modes to consider (in rows of m, n). Used to generate considered TM modes.
        
        interpolationPointsTau = 2^12;  % Number of points to use for the interpolation function along tau.
        integralPointsTauLossy = 135;   % Number of points to integrate along tau for lossy structures.
        integralPointsPsi = 50;         % Number of points to use to integrate along psi.

        convergenceAbsTol = 0.001;      % Default tolerance for reflection coefficient calculations.
    end
    
    %% Public member function definitions (implemented in separate files)
    methods (Access = public)
        gam = calculate(O, f, er, ur, thk, AbsTol);
        
        [modes] = setModes(O, maxM, maxN);
        [a, b] = setWaveguideBand(O, band);
        setWaveguideDimensions(O, a, b);
        recomputeInterpolants(O);
    end
    
    %% Private member function definitions (implemented in separate files)
    methods (Access = private)
        [A1b1] = integralFun(O, tauP, er, ur, thk, k0);
        [tau, weightsE, weightsH] = computeIntegralInterpolant(O, orderTau);
        [A1, A2, b1, b2] = constructMatrixEquation(O, nLayerInt);
        [k_A1, k_A2, k_b1, k_b2] = constructFrequencyMultipliers(O, k0);
    end
    
    %% Private static function definitions (implemented in separate files)
    methods (Static, Access = public)
        [specE, specH] = multilayerSpectrumRect(tau, k0, er, ur, thk);
    end
    
    %% Class constructor
    methods
        function O = nLayerRectangular(varargin)
            %NLAYERRECTANGULAR Construct an instance of this class
            %   Usage:
            %       nLayerRectangular(a, b, maxM, maxN);
            %       nLayerRectangular(band, maxM, maxN);
            %   
            %   When calling as nLayerRectangular(a, b, maxM, maxN), the 
            %   inputs a and b are the broad and narrow dimensions,
            %   respectively, in mm, and maxM and maxN are the highest
            %   mode indices to consider.
            %
            %   When calling as nLayerRectangular(band, maxM, maxN), band
            %   is a case-insensitive waveguide band identifier, and and 
            %   maxM and maxN are the highest mode indices to consider.
            
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

