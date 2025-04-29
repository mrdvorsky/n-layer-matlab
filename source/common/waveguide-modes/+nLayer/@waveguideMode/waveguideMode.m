classdef waveguideMode < matlab.mixin.Copyable
    %Class defining properties of a waveguide mode.
    % This class defines the properties of a particular mode in a
    % waveguide, including the mode type, symmetries, spectrums, etc.
    %
    % Author: Matt Dvorsky

    properties (Access=public)
        modeType {mustBeTextScalar, mustBeMember(modeType, ["TE", "TM"])} = "TE";   % Type of mode (TE or TM).
        modeLabel {mustBeTextScalar} = "";      % Human readable mode name.

        kc0(1, 1) {mustBeNonnegative} = 0;      % Mode cutoff wavenumber when filled with vacuum.
        apertureSize(1, 1) {mustBePositive} = 1;   % Furthest distance between two points on aperture.

        symmetryX {mustBeTextScalar, mustBeMember(...   % Symmetry of mode about x-axis.
            symmetryX, ["None", "PEC", "PMC"])} = "None";
        symmetryY {mustBeTextScalar, mustBeMember(...   % Symmetry of mode about y-axis.
            symmetryY, ["None", "PEC", "PMC"])} = "None";
        symmetryAxial {mustBeTextScalar, mustBeMember(...   % Axial symmetry about origin.
            symmetryAxial, ["None", "TE", "TM"])} = "None";

        WhSpec(1, 1) ...    % Fourier spectrum function for TE portion of waveguide mode.
            {mustBeCallable(WhSpec, {1, 0, 1, 0}, "kx, ky, kr, kphi")} = @(~, ~, ~, ~) 0;
        WeSpec(1, 1) ...    % Fourier spectrum function for TM portion of waveguide mode.
            {mustBeCallable(WeSpec, {1, 0, 1, 0}, "kx, ky, kr, kphi")} = @(~, ~, ~, ~) 0;
    end
    properties(Dependent, GetAccess=public, SetAccess=protected)
        ExSpec(1, 1);   % Fourier spectrum function for x-component of electric field.
        EySpec(1, 1);   % Fourier spectrum function for y-component of electric field.
        EzSpec(1, 1);   % Fourier spectrum function for z-component of electric field.
    end

    %% Class Constructor
    methods
        function O = waveguideMode(classProperties)
            %Construct an instance of this class.

            arguments
                classProperties.?nLayer.waveguideMode;
            end

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                O.(propPairs{ii}) = propPairs{ii + 1};
            end
        end
    end

    %% Class Functions
    methods (Access=public)
        [] = showMode(O, options);
    end

    %% Class Getters
    methods
        function [spec] = get.ExSpec(O)
            if strcmp(O.modeType, "TE")
                spec = @(kx, ky, kr, kphi) ...
                    cos(kphi).*O.WeSpec(kx, ky, kr, kphi) ...
                    + sin(kphi).*O.WhSpec(kx, ky, kr, kphi);
            else
                spec = @(kx, ky, kr, kphi) ...
                    cos(kphi).*O.WeSpec(kx, ky, kr, kphi);
            end
        end
        function [spec] = get.EySpec(O)
            if strcmp(O.modeType, "TE")
                spec = @(kx, ky, kr, kphi) ...
                    sin(kphi).*O.WeSpec(kx, ky, kr, kphi) ...
                    - cos(kphi).*O.WhSpec(kx, ky, kr, kphi);
            else
                spec = @(kx, ky, kr, kphi) ...
                    cos(kphi).*O.WeSpec(kx, ky, kr, kphi);
            end
        end
        function [spec] = get.EzSpec(O)
            if strcmp(O.modeType, "TE")
                spec = @(~, ~, ~, ~) 0;
            else
                spec = @(kx, ky, kr, kphi) ...
                    O.WeSpec(kx, ky, kr, kphi) ./ kr;
            end
        end
    end
end

