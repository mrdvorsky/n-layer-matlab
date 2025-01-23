classdef rectangularWaveguideBand
    %RECTANGULARWAVEGUIDEBAND Defines waveguide bands and their dimensions.
    % Contains the definitions of standard rectangular waveguide bands. To
    % add a new band, edit the enumeration section below.
    %
    % The intended usage of this class is either as an argument validator
    % for a function that accepts a waveguide band, or simply as a function
    % to get waveguide dimensions. See example usage below.
    %
    % Example Usage:
    %   arguments
    %       ...
    %       wgBand(1, 1) rectangularWaveguideBand;  % Accepts one band.
    %       wgBands(:, 1) rectangularWaveguideBand; % Accepts an array of bands.
    %       ...
    %   end
    %
    %   wgBand = rectangularWaveguideBand.Ka;       % Get Ka-band object.
    %   [a, b] = wgBand.getDimensions();            % Get dimensions in mm.
    %   [a, b] = wgBand.getDimensions(0.001);       % Get dimensions in mm.
    %   [a, b] = wgBand.getDimensions(1);           % Get dimensions in m.
    %   [a, b] = wgBand.getDimensions(scale);       % Get dimensions in m / scale.
    %
    % Author: Matt Dvorsky

    enumeration     % Define as "BandName  (a_inches, b_inches)".
        R      (4.300, 2.150);
        WR430  (4.300, 2.150);
        S      (2.840, 1.340);
        WR284  (2.840, 1.340);
        G      (1.872, 0.872);
        WR187  (1.872, 0.872);
        J      (1.372, 0.622);
        WR137  (1.372, 0.622);
        X      (0.900, 0.400);
        WR90   (0.900, 0.400);
        Ku     (0.622, 0.311);
        WR62   (0.622, 0.311);
        K      (0.420, 0.170);
        WR42   (0.420, 0.170);
        Ka     (0.280, 0.140);
        WR28   (0.280, 0.140);
        Q      (0.224, 0.112);
        WR22   (0.224, 0.112);
        U      (0.188, 0.094);
        WR19   (0.188, 0.094);
        V      (0.148, 0.074);
        WR15   (0.148, 0.074);
        E      (0.122, 0.061);
        WR12   (0.122, 0.061);
        W      (0.100, 0.050);
        WR10   (0.100, 0.050);

        Custom (0, 0);
    end

    %% Class Implementation
    properties
        waveguideA;
        waveguideB;
    end
    methods (Access=public)
        function O = rectangularWaveguideBand(a_inches, b_inches)
            O.waveguideA = 0.0254 * a_inches;
            O.waveguideB = 0.0254 * b_inches;
        end
    end
    methods (Access=public)
        function [a, b] = getDimensions(O, distanceScaleFactor)
            arguments
                O;
                distanceScaleFactor(1, 1) {mustBePositive, mustBeFinite} = 0.001;
            end
            a = O.waveguideA ./ distanceScaleFactor;
            b = O.waveguideB ./ distanceScaleFactor;
        end
    end
end

