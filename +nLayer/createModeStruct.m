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
    modeStruct.CutoffWavenumber(1, 1) {mustBeNonnegative};
    modeStruct.MaxOperatingWavenumber(1, 1) {mustBeNonnegative};

    modeStruct.IsExcitationMode(1, 1) logical = false;
    modeStruct.IsReceiveMode(1, 1) logical = false;

    modeStruct.SymmetryX {mustBeTextScalar, mustBeMember(...
        modeStruct.SymmetryX, ["None", "Even", "Odd"])} = "None";
    modeStruct.SymmetryY {mustBeTextScalar, mustBeMember(...
        modeStruct.SymmetryY, ["None", "Even", "Odd"])} = "None";
    modeStruct.SymmetryAxial {mustBeTextScalar, mustBeMember(...
        modeStruct.SymmetryAxial, ["None", "TE", "TM"])} = "None";

    modeStruct.ApertureWidth(1, 1) {mustBePositive} = 1;

    modeStruct.WaveguideEr(1, 1) {nLayer.mustBeErUrCallable} = @(f) 1;
    modeStruct.WaveguideUr(1, 1) {nLayer.mustBeErUrCallable} = @(f) 1;

    modeStruct.OffsetX(1, 1) {mustBeReal, mustBeFinite} = 0;
    modeStruct.OffsetY(1, 1) {mustBeReal, mustBeFinite} = 0;
    modeStruct.RotationAngle(1, 1) {mustBeReal, mustBeFinite} = 0;
end

%% Assign Mode Type and Label
modeStruct.ModeType = modeType;
modeStruct.ModeLabel = modeLabel;

%% Check Mode Spectrums
if ~isfield(modeStruct, "ExSpec") && ~isfield(modeStruct, "EySpec")
    error("At least one electric field spectrum (i.e., 'ExSpec', " + ...
        "'EySpec') must be specified.");
end

if ~isfield(modeStruct, "ExSpec")
    modeStruct.ExSpec = @(~, ~, ~, ~) 0;
end

if ~isfield(modeStruct, "EySpec")
    modeStruct.EySpec = @(~, ~, ~, ~) 0;
end

%% Check WaveguideEr and Ur
if isnumeric(modeStruct.WaveguideEr)
    er = modeStruct.WaveguideEr;
    modeStruct.WaveguideEr = @(f) er;
end

if isnumeric(modeStruct.WaveguideUr)
    ur = modeStruct.WaveguideUr;
    modeStruct.WaveguideUr = @(f) ur;
end

%% Check Mode Cutoff Wavenumber
if ~isfield(modeStruct, "CutoffWavenumber")
    error("Cutoff wavenumber (i.e., 'CutoffWavenumber') must be specified.");
end

if ~isfield(modeStruct, "MaxOperatingWavenumber")
    modeStruct.MaxOperatingWavenumber = 2*modeStruct.CutoffWavenumber;
    warning("Maximum operating frequency (i.e., 'MaxOperatingWavenumber') " + ...
        "was not specified for the mode. A default value of " + ...
        "'2 * CutoffWavenumber' will be used.");
end

%% Check Symmetries
if strcmp(modeStruct.SymmetryAxial, "TE")
    modeStruct.SymmetryX = "Odd";
    modeStruct.SymmetryY = "Odd";
elseif strcmp(modeStruct.SymmetryAxial, "TM")
    modeStruct.SymmetryX = "Even";
    modeStruct.SymmetryY = "Even";
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
    "CutoffWavenumber", ...
    "MaxOperatingWavenumber", ...
    "ApertureWidth", ...
    "WaveguideEr", ...
    "WaveguideUr", ...
    "OffsetX", ...
    "OffsetY", ...
    "RotationAngle"]);

end




