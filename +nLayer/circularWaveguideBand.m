classdef circularWaveguideBand
    %CIRCULARWAVEGUIDEBAND Defines waveguide bands and their dimensions.
    % Contains the definitions of standard circular waveguide bands, in
    % addition to the custom TE01 waveguide bands. To add a new band, edit
    % the enumeration section below.
    %
    % The intended usage of this class is either as an argument validator
    % for a function that accepts a waveguide band, or simply as a function
    % to get waveguide dimensions. See example usage below.
    %
    % Example Usage:
    %   arguments
    %       ...
    %       wgBand(1, 1) circularWaveguideBand;  % Accepts one band.
    %       wgBands(:, 1) circularWaveguideBand; % Accepts an array of bands.
    %       ...
    %   end
    %
    %   wgBand = circularWaveguideBand.KaTE01;      % Get Ka-band object.
    %   [r] = wgBand.getDimensions();               % Get radius in mm.
    %   [r] = wgBand.getDimensions(0.001);          % Get radius in mm.
    %   [r] = wgBand.getDimensions(1);              % Get radius in m.
    %   [r] = wgBand.getDimensions(scale);          % Get radius in m / scale.
    %
    % Author: Matt Dvorsky

    enumeration     % Define as "BandName  (r_inches)".
        X_TE01   (93/128);
        Ka_TE01  (30/128);
        Q_TE01   (23/128);
        V_TE01   (17/128);

        X_Low    (1.094);
        X_Mid    (0.938);
        X_High   (0.797);
        Ku_Low   (0.688);
        Ku_Mid   (0.594);
        Ku_High  (0.500);
        K_Low    (0.455);
        K_Mid    (0.396);
        K_High   (0.328);
        Ka_Low   (0.315);
        Ka_Mid   (0.250);
        Ka_High  (0.219);
        Q_Low    (0.250);
        Q_Mid    (0.219);
        Q_High   (0.188);
        U_Low    (0.210);
        U_Mid    (0.188);
        U_High   (0.165);
        V_Low    (0.165);
        V_Mid    (0.141);
        V_High   (0.125);
        E_Low    (0.136);
        E_Mid    (0.125);
        E_High   (0.094);
        W_Low    (0.112);
        W_High   (0.094);
        F_Low    (0.089);
        F_High   (0.075);
        D_Low    (0.073);
        D_High   (0.059);
        G_Low    (0.058);
        G_High   (0.045);

        Custom   (0);
    end

    %% Class Implementation
    properties
        waveguideR;
    end
    methods (Access=public)
        function O = circularWaveguideBand(r_inches)
            O.waveguideR = 0.0254 * r_inches;
        end
    end
    methods (Access=public)
        function [r] = getDimensions(O, distanceScaleFactor)
            arguments
                O;
                distanceScaleFactor(1, 1) {mustBePositive, mustBeFinite} = 0.001;
            end
            r = O.waveguideR ./ distanceScaleFactor;
        end
    end
end

