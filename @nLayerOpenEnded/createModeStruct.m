function [modeStruct] = createModeStruct(modeStruct)
%CREATEMODESTRUCT Creates and formats a mode structure.
% This function creates a mode structure that defines the modes of a
% waveguide. Pass the mode spectrums, cutoff betas, symmetry conditions,
% etc., to this function to create the structure.
%
% This function should be called in the "defineWaveguideModes" function
% that is overloaded by a subclass.
%
% Outputs:
%   modeStruct - Mode struct with information about waveguide modes.
%
% Named Arguments:
%   PositiveOddInt - Positive odd integer optional named argument. No default.
%   StringArray ("") - String array named optional argument. Default
%       value is an array with an empty string and should be listed in
%       parentheses.
%
% Author: Matt Dvorsky

arguments
    modeStruct.SpecEx_TE(:, 1) {mustBeA(modeStruct.SpecEx_TE, "cell")} = {};
    modeStruct.SpecEy_TE(:, 1) {mustBeA(modeStruct.SpecEy_TE, "cell")} = {};
    modeStruct.SpecEx_TM(:, 1) {mustBeA(modeStruct.SpecEx_TM, "cell")} = {};
    modeStruct.SpecEy_TM(:, 1) {mustBeA(modeStruct.SpecEy_TM, "cell")} = {};
    modeStruct.SpecEx_Hybrid(:, 1) {mustBeA(modeStruct.SpecEx_Hybrid, "cell")} = {};
    modeStruct.SpecEy_Hybrid(:, 1) {mustBeA(modeStruct.SpecEy_Hybrid, "cell")} = {};

    modeStruct.CutoffBeta_TE(:, 1) = [];
    modeStruct.CutoffBeta_TM(:, 1) = [];
    modeStruct.CutoffBeta_Hybrid(:, 1) = [];

    modeStruct.ModeSymmetryX {mustBeMember(modeStruct.ModeSymmetryX, ...
        ["None", "Even", "Odd"])} = "None";
    modeStruct.ModeSymmetryY {mustBeMember(modeStruct.ModeSymmetryY, ...
        ["None", "Even", "Odd"])} = "None";
    modeStruct.ModeSymmetryAxial {mustBeMember(modeStruct.ModeSymmetryAxial, ...
        ["None", "TE", "TM"])} = "None";

    modeStruct.IntegralScaleFactor(1, 1) {mustBePositive} = 1;
end

%% Check Inputs
% Check sizes of mode spectrums.
[modeStruct.SpecEx_TE, modeStruct.SpecEy_TE] = ...
    checkSpectrumSizes(modeStruct.SpecEx_TE, modeStruct.SpecEy_TE);
[modeStruct.SpecEx_TM, modeStruct.SpecEy_TM] = ...
    checkSpectrumSizes(modeStruct.SpecEx_TM, modeStruct.SpecEy_TM);
[modeStruct.SpecEx_Hybrid, modeStruct.SpecEy_Hybrid] = ...
    checkSpectrumSizes(modeStruct.SpecEx_Hybrid, modeStruct.SpecEy_Hybrid);

% Check cutoff beta definitions.
if (numel(modeStruct.CutoffBeta_TE) ~= numel(modeStruct.SpecEx_TE)) ...
        || (numel(modeStruct.CutoffBeta_TM) ~= numel(modeStruct.SpecEx_TM)) ...
        || (numel(modeStruct.CutoffBeta_Hybrid) ~= numel(modeStruct.SpecEx_Hybrid))
    error("CutoffBeta arguments size must match size of " + ...
        "corresponding mode spectrums");
end


end


%% Input Argument Checking
function [specEx, specEy] = checkSpectrumSizes(specEx, specEy)

if isempty(specEx) && ~isempty(specEy)
    specEx = repmat({@(~, ~, ~, ~) 0}, numel(specEy), 1);
end
if isempty(specEy) && ~isempty(specEx)
    specEy = repmat({@(~, ~, ~, ~) 0}, numel(specEx), 1);
end
if numel(specEx) ~= numel(specEy)
    error("Mode spectrum arguments for Ex and Ey must " + ...
        "be cell arrays of the same size or one must be empty.");
end

end

