function [modeStruct] = createModeStruct(modeType, modeLabel, modeStruct)
%CREATEMODESTRUCT Creates and formats a mode structure for a waveguide mode.
% This function creates a mode structure that defines a modes in a
% waveguide. Pass the in the electric field spectrum functions, cutoff
% wavenumbers, symmetry conditions, etc., to this function to create the
% structure.
%
% This function should be called in the "defineWaveguideModes" function
% that is overloaded by a subclass.
%
% Outputs:
%   modeStruct - Mode struct with information about the waveguide mode.
%
% Named Arguments:
%   PositiveOddInt - .
%
% Author: Matt Dvorsky

arguments
    modeType {mustBeTextScalar, mustBeMember(modeType, ...
        ["TE", "TM", "Hybrid"])};
    modeLabel {mustBeTextScalar};

    modeStruct.ExSpec(1, 1) {mustBeCallable(...
        modeStruct.ExSpec, {1, 0, 1, 0}, "kx, ky, kr, kphi")};
    modeStruct.EySpec(1, 1) {...
        mustBeCallable(modeStruct.EySpec, {1, 0, 1, 0}, "kx, ky, kr, kphi")};

    modeStruct.WhSpec(1, 1) {mustBeCallable(...
        modeStruct.WhSpec, {1, 0, 1, 0}, "kx, ky, kr, kphi")};
    modeStruct.WeSpec(1, 1) {...
        mustBeCallable(modeStruct.WeSpec, {1, 0, 1, 0}, "kx, ky, kr, kphi")};

    modeStruct.CutoffWavenumber(1, 1) {mustBeNonnegative};

    modeStruct.IsExcitationMode(1, 1) logical = false;
    modeStruct.IsReceiveMode(1, 1) logical = false;

    modeStruct.SymmetryX string {mustBeMember(...
        modeStruct.SymmetryX, ["PEC", "PMC", "None"])} = "None";
    modeStruct.SymmetryY string {mustBeMember(...
        modeStruct.SymmetryY, ["PEC", "PMC", "None"])} = "None";
    modeStruct.SymmetryAxial string {mustBeMember(...
        modeStruct.SymmetryAxial, ["TE", "TM", "None"])} = "None";

    modeStruct.ApertureWidth(1, 1) {mustBePositive} = 1;

    modeStruct.WaveguideEr(1, 1) {nLayer.mustBeErUrCallable} = 1;
    modeStruct.WaveguideUr(1, 1) {nLayer.mustBeErUrCallable} = 1;

    modeStruct.OffsetX(1, 1) {mustBeReal, mustBeFinite} = 0;
    modeStruct.OffsetY(1, 1) {mustBeReal, mustBeFinite} = 0;
    modeStruct.RotationAngle(1, 1) {mustBeReal, mustBeFinite} = 0;

    modeStruct.HertzAz(1, 1) {...
        mustBeCallable(modeStruct.HertzAz, {1, 1}, "x, y")} = @(x, y) 0;
    modeStruct.HertzAz_dx(1, 1) {...
        mustBeCallable(modeStruct.HertzAz_dx, {1, 1}, "x, y")} = @(x, y) 0;
    modeStruct.HertzAz_dy(1, 1) {...
        mustBeCallable(modeStruct.HertzAz_dy, {1, 1}, "x, y")} = @(x, y) 0;

    modeStruct.BoundaryPoints(:, 2) {mustBeReal, mustBeFinite} = [];
end

%% Assign Mode Type and Label
modeStruct.ModeType = modeType;
modeStruct.ModeLabel = modeLabel;

%% Check Mode Spectrums
kc0 = modeStruct.CutoffWavenumber;

WhSpec = modeStruct.WhSpec;
WeSpec = modeStruct.WeSpec;
if strcmp(modeStruct.ModeType, "TE")
    modeStruct.ExSpec = @(kx, ky, kr, kphi) ...
        cos(kphi).*WeSpec(kx, ky, kr, kphi) ...
        + sin(kphi).*WhSpec(kx, ky, kr, kphi);
    modeStruct.EySpec = @(kx, ky, kr, kphi) ...
        sin(kphi).*WeSpec(kx, ky, kr, kphi) ...
        - cos(kphi).*WhSpec(kx, ky, kr, kphi);
    modeStruct.EzSpec = @(kx, ky, kr, kphi) ...
        WhSpec(kx, ky, kr, kphi) .* (kr ./ kc0.^2);
else
    modeStruct.ExSpec = @(kx, ky, kr, kphi) ...
        cos(kphi).*WeSpec(kx, ky, kr, kphi);
    modeStruct.EySpec = @(kx, ky, kr, kphi) ...
        sin(kphi).*WeSpec(kx, ky, kr, kphi);
    modeStruct.EzSpec = @(kx, ky, kr, kphi) ...
        WeSpec(kx, ky, kr, kphi) ./ kr;
end


%% Check Mode Cutoff Wavenumber
if ~isfield(modeStruct, "CutoffWavenumber")
    error("Cutoff wavenumber (i.e., 'CutoffWavenumber') must be specified.");
end

if ~isfield(modeStruct, "MaxOperatingWavenumber")
    modeStruct.MaxOperatingWavenumber = 2*modeStruct.CutoffWavenumber;
    % warning("Maximum operating frequency (i.e., 'MaxOperatingWavenumber') " + ...
    %     "was not specified for the mode. A default value of " + ...
    %     "'2 * CutoffWavenumber' will be used.");
end

%% Check Symmetries
if strcmp(modeStruct.SymmetryAxial, "TE")
    modeStruct.SymmetryX = "PEC";
    modeStruct.SymmetryY = "PEC";
elseif strcmp(modeStruct.SymmetryAxial, "TM")
    modeStruct.SymmetryX = "PMC";
    modeStruct.SymmetryY = "PMC";
end

%% Check Offsets and Rotations
if modeStruct.OffsetX ~= 0
    warning("'OffsetX' is currently unsupported.");
end
if modeStruct.OffsetY ~= 0
    warning("'OffsetY' is currently unsupported.");
end
if modeStruct.RotationAngle ~= 0
    warning("'RotationAngle' is currently unsupported.");
end

%% Order Field Names
modeStruct = orderfields(modeStruct, [...
    "ModeType", ...
    "ModeLabel", ...
    "IsExcitationMode", ...
    "IsReceiveMode", ...
    "SymmetryX", ...
    "SymmetryY", ...
    "SymmetryAxial", ...
    "ExSpec", ...
    "EySpec", ...
    "EzSpec", ...
    "WhSpec", ...
    "WeSpec", ...
    "CutoffWavenumber", ...
    "MaxOperatingWavenumber", ...
    "ApertureWidth", ...
    "WaveguideEr", ...
    "WaveguideUr", ...
    "OffsetX", ...
    "OffsetY", ...
    "RotationAngle", ...
    "HertzAz", ...
    "HertzAz_dx", ...
    "HertzAz_dy", ...
    "BoundaryPoints"]);

end




