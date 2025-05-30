classdef coaxialConnectorType
    %Defines coaxial connector types and their dimensions.
    % Contains the definitions of standard coaxial waveguide connectors.
    % To add a new connector type, edit the enumeration section below.
    %
    % The intended usage of this class is either as an argument validator
    % for a function that accepts a coaxial connector type, or simply as a
    % function to get coax dimensions. See example usage below.
    %
    % Example Usage:
    %   arguments
    %       ...
    %       coaxType(1, 1)  nLayer.coaxialConnectorType;   % Accepts one connector.
    %       coaxTypes(:, 1) nLayer.coaxialConnectorType;   % Accepts an array of connectors.
    %       ...
    %   end
    %
    %   coaxType = coaxialConnectorType.N;          % Get N-type object.
    %   [di, do] = coaxType.getDimensions();        % Get dimensions in mm.
    %   [di, do] = coaxType.getDimensions(0.001);   % Get dimensions in mm.
    %   [di, do] = coaxType.getDimensions(1);       % Get dimensions in m.
    %   [di, do] = coaxType.getDimensions(scale);   % Get dimensions in m / scale.
    %   er = coaxType.fillEr;                       % Get coax fill permittivity (er).
    %
    % Author: Matt Dvorsky

    enumeration     % Define as "CoaxTypeName  (do_mm, di_mm, er)".
        N      (7.000, 3.039, 1);
        K      (2.920, 1.268, 1);

        Custom (1, 0, 1);
    end

    %% Class Implementation
    properties
        innerDiameter;
        outerDiameter;
        fillEr;
    end
    methods (Access=public)
        function self = coaxialConnectorType(do_mm, di_mm, er)
            if di_mm >= do_mm
                error("Outer radius must be larger than inner.");
            end
            self.innerDiameter = 0.001 * di_mm;
            self.outerDiameter = 0.001 * do_mm;
            self.fillEr = er;
        end
    end
    methods (Access=public)
        function [di, do] = getDimensions(self, distanceScaleFactor)
            arguments
                self nLayer.coaxialConnectorType;
                distanceScaleFactor(1, 1) {mustBePositive, mustBeFinite} = 0.001;
            end
            di = self.innerDiameter ./ distanceScaleFactor;
            do = self.outerDiameter ./ distanceScaleFactor;
        end
    end
end

