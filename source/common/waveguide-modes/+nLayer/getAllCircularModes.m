function [waveguideModes] = getAllCircularModes(m, n, wgR, options)
%Get all possible "waveguideMode" objects for a circular waveguide.
% This functions returns "waveguideMode" objects for all circular
% waveguide modes that match the pattern TEmn or TMmn, for all
% combinations of the input vectors "m" and "n";
%
% Optionally, symmetry filters can be applied, so that all returned mode
% objects have the specified symmetry.
%
% Example Usage:
%   % All modes, regardless of symmetry.
%   [waveguideModes] = getAllCircularModes(m, n, wgR);
%
%   % Only return modes where the x-axis could be replaced with PEC.
%   [waveguideModes] = getAllCircularModes(m, n, wgR, modeSymmetryX="PEC");
%
%
% Inputs:
%   m - Vector of "m" values for returned TEmn and TMmn modes.
%   n - Vector of "n" values for returned TEmn and TMmn modes.
%   wgR - Radius of circular waveguide.
%
% Author: Matt Dvorsky

arguments
    m(:, 1) {mustBeInteger, mustBeNonnegative};
    n(1, :) {mustBeInteger, mustBePositive};
    wgR(1, 1) {mustBePositive};

    options.SymmetryX string {mustBeMember(options.SymmetryX, ...
        ["PEC", "PMC", "None"])} = "None";
    options.SymmetryY string {mustBeMember(options.SymmetryY, ...
        ["PEC", "PMC", "None"])} = "None";
    options.SymmetryAxial string {mustBeMember(options.SymmetryAxial, ...
        ["TE", "TM", "None"])} = "None";
end

%% Generate List of Modes
if strcmp(options.SymmetryAxial, "TE")
    modes_TE = [0*(1:n); (1:n)].';
    modes_TM = [];
elseif strcmp(options.SymmetryAxial, "TM")
    modes_TM = [0*(1:n); (1:n)].';
    modes_TE = [];
else
    modes_TE = [reshape((m).' + 0*(1:n), [], 1), ...
        reshape(0*(m).' + (1:n), [], 1)];
    modes_TM = modes_TE;
end

%% Get "nLayer.waveguideMode" Objects
modesAll = [modes_TE; modes_TM];
modeTypes = [repmat("TE", size(modes_TE, 1), 1); ...
    repmat("TM", size(modes_TM, 1), 1)];

isRotated = false(size(modesAll, 1), 1);
modeTypes = [modeTypes; modeTypes(modesAll(:, 1) ~= 0)];
modesAll = [modesAll; modesAll(modesAll(:, 1) ~= 0, :)];
isRotated = [isRotated; true(size(modesAll, 1) - numel(isRotated), 1)];

%#ok<*AGROW>
waveguideModes = nLayer.waveguideMode.empty;
for ii = 1:size(modesAll, 1)
    m = modesAll(ii, 1);
    n = modesAll(ii, 2);
    waveguideModes(1, ii) = nLayer.getCircularModeStruct(...
        m, n, wgR, modeTypes(ii), isRotated(ii));
end

%% Sort by Cutoff
[~, sortInd] = sort([waveguideModes.kc0]);
waveguideModes = waveguideModes(sortInd);

end

