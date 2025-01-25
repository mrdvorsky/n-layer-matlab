function [waveguideModes] = getCircularModes(m, maxN, wgR, options)
%GETCIRCULARMODES Get "waveguideMode" objects for a circular waveguide.
%
% Author: Matt Dvorsky

arguments
    m(1, 1) {mustBeInteger, mustBeNonnegative};
    maxN(1, 1) {mustBeInteger, mustBePositive};
    wgR(1, 1) {mustBePositive};

    options.SymmetryX string {mustBeMember(options.SymmetryX, ...
        ["PEC", "PMC", "None"])} = "PEC";
    options.SymmetryY string {mustBeMember(options.SymmetryY, ...
        ["PEC", "PMC", "None"])} = "PMC";
    options.SymmetryAxial string {mustBeMember(options.SymmetryAxial, ...
        ["TE", "TM", "None"])} = "None";
end

%% Generate List of Modes
if strcmp(options.SymmetryAxial, "TE")
    modes_TE = [0*(1:maxN); (1:maxN)].';
    modes_TM = [];
elseif strcmp(options.SymmetryAxial, "TM")
    modes_TM = [0*(1:maxN); (1:maxN)].';
    modes_TE = [];
else
    modes_TE = [reshape((m).' + 0*(1:maxN), [], 1), ...
        reshape(0*(m).' + (1:maxN), [], 1)];
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

