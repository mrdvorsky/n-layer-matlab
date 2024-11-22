classdef waveguideMode
    %WAVEGUIDEMODE Class defining properties of a waveguide mode.
    % This class defines the properties of a particular mode in a
    % waveguide, including the mode type, symmetries, spectrums, etc.
    %
    % Author: Matt Dvorsky

    properties
        modeType {mustBeTextScalar, mustBeMember(modeType, ["TE", "TM", "Hybrid"])};
        modeLabel {mustBeTextScalar};

        ExSpec(1, 1) {mustBeCallable(ExSpec, {1, 0, 1, 0}, "kx, ky, kr, kphi")};
        EySpec(1, 1) {mustBeCallable(EySpec, {1, 0, 1, 0}, "kx, ky, kr, kphi")};
        CutoffWavenumber(1, 1) {mustBeNonnegative};
        MaxOperatingWavenumber(1, 1) {mustBeNonnegative};

        IsExcitationMode(1, 1) logical = false;
        IsReceiveMode(1, 1) logical = false;

        SymmetryX {mustBeTextScalar, mustBeMember(...
            SymmetryX, ["None", "Even", "Odd"])} = "None";
        SymmetryY {mustBeTextScalar, mustBeMember(...
            SymmetryY, ["None", "Even", "Odd"])} = "None";
        SymmetryAxial {mustBeTextScalar, mustBeMember(...
            SymmetryAxial, ["None", "TE", "TM"])} = "None";

        apertureWidth(1, 1) {mustBePositive} = 1;

        waveguideEr(1, 1) {nLayer.mustBeErUrCallable} = @(f) 1;
        waveguideUr(1, 1) {nLayer.mustBeErUrCallable} = @(f) 1;

        offsetX(1, 1) {mustBeReal, mustBeFinite} = 0;
        offsetY(1, 1) {mustBeReal, mustBeFinite} = 0;
        rotationAngle(1, 1) {mustBeReal, mustBeFinite} = 0;
    end

    methods
        function obj = waveguideMode(inputArg1,inputArg2)
            %WAVEGUIDEMODE Construct an instance of this class.
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

