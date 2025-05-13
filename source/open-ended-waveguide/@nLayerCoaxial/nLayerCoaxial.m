classdef nLayerCoaxial < nLayerOpenEnded
    %Implementation of nLayerForward for coaxial waveguides.
    % This class can be used to calculate the reflection coefficient seen
    % by a coaxial waveguide looking into a multilayer structure. Note
    % that the units of all parameters should match that of the speed of
    % light specified by the "distanceUnitScale" and "frequencyUnitScale"
    % parameters (defaults are mm and GHz), both of which can be set in
    % the constructor.
    %
    % Example Usage:
    %   NL = nLayerCoaxial(m, maxN, waveguideBand="N");
    %   NL = nLayerCoaxial(m, maxN, waveguideRi=0.44, waveguideRi=1);
    %   NL = nLayerCoaxial(m, maxN, distanceUnitScale=1, ...
    %       waveguideRi=0.44e-3, waveguideRi=1e-3);
    %   NL = nLayerCoaxial(m, maxN, waveguideBand="K", ...
    %       prop1=val1, prop2=val2, ...);
    %
    %   gam = NL.calculate(f, er, ur, thk);
    %   gam = NL.calculate(f, er, {}, thk);
    %   gam = NL.calculate(f, {}, ur, thk);
    %
    %
    % Author: Matt Dvorsky

    properties (Access=public, AbortSet)
        waveguideBand(1, 1) nLayer.coaxialConnectorType;    % Waveguide band.

        modeIndexM(1, 1) {mustBeInteger, mustBeNonnegative} = 0;    % Value of 'm' for considered TEmn or TMmn modes.
        maxModeIndexN(1, 1) {mustBeInteger, mustBeNonnegative} = 0; % Maximum value of 'n' for considered TEmn or TMmn modes.
    end
    properties (Dependent, Access=public, AbortSet)
        waveguideRi(1, 1) {mustBePositive, mustBeFinite};   % Coax inner radius.
        waveguideRo(1, 1) {mustBePositive, mustBeFinite};   % Coax outer radius.
    end
    properties (Access=private)
        waveguideRi_custom(1, 1);
        waveguideRo_custom(1, 1);
    end

    %% Class Constructor
    methods
        function self = nLayerCoaxial(indexM, maxIndexN, classProperties)
            %Construct an instance of this class.
            % Inputs:
            %   indexM - Index 'm' of TEmn and TMmn modes to consider.
            %   maxIndexN - Highest index 'n' of TEmn and TMmn modes to consider.

            arguments
                indexM(1, 1) {mustBeInteger, mustBeNonnegative};
                maxIndexN(1, 1) {mustBeInteger, mustBeNonnegative};
                classProperties.?nLayerCoaxial;
            end

            self.modeIndexM = indexM;
            self.maxModeIndexN = maxIndexN;

            if self.modeIndexM == 0
                self.modeSymmetryAxial = "TM";
            end

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                self.(propPairs{ii}) = propPairs{ii + 1};
            end

            % if strcmp(O.modeSymmetryAxial, "TM")
            %     [kc] = besselCrossZeros(O.modeIndexM, O.waveguideRo./O.waveguideRi, 1);
            %     O.frequencyRange = [0, kc] ./ (2*pi) ...
            %         .* O.speedOfLight;
            % end

            if isfield(classProperties, "frequencyRange")
                self.frequencyRange = classProperties.frequencyRange;
            end
        end
    end

    %% Class Functions
    methods (Access=protected)
        [waveguideModes] = defineWaveguideModes(self);
    end

    %% Class Setters
    methods
        function set.waveguideBand(self, newBand)
            self.waveguideBand = newBand;
            self.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.waveguideRi(self, newA)
            self.waveguideRi_custom = newA;
            self.waveguideRo_custom = self.waveguideRo;
            self.waveguideBand = "Custom";
        end
        function set.waveguideRo(self, newB)
            self.waveguideRo_custom = newB;
            self.waveguideRi_custom = self.waveguideRi;
            self.waveguideBand = "Custom";
        end

        function set.modeIndexM(self, newMaxInd)
            self.modeIndexM = newMaxInd;
            self.shouldRegenerateWaveguideModeObjects = true;
        end
        function set.maxModeIndexN(self, newMaxInd)
            self.maxModeIndexN = newMaxInd;
            self.shouldRegenerateWaveguideModeObjects = true;
        end
    end

    %% Class Getters
    methods
        function [ri] = get.waveguideRi(self)
            if self.waveguideBand == "Custom"
                ri = self.waveguideRi_custom;
                return;
            end
            [di, ~] = self.waveguideBand.getDimensions(self.distanceUnitScale);
            ri = 0.5 * di;
        end
        function [ro] = get.waveguideRo(self)
            if self.waveguideBand == "Custom"
                ro = self.waveguideRo_custom;
                return;
            end
            [~, do] = self.waveguideBand.getDimensions(self.distanceUnitScale);
            ro = 0.5 * do;
        end
    end

end

