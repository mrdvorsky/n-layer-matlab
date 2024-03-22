function [modeStructs] = getCoaxialModes(maxM, maxN, wgRi, wgRo, options)
%GETCOAXIALMODES Get "waveguideMode" objects for coaxial waveguide.
%
% Author: Matt Dvorsky

arguments
    maxM(1, 1) {mustBeInteger, mustBeNonnegative};
    maxN(1, 1) {mustBeInteger, mustBeNonnegative};
    wgRi(1, 1) {mustBePositive};
    wgRo(1, 1) {mustBePositive};

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
    modes_TM = [0*(0:maxN); (0:maxN)].';
    modes_TE = [];
else
    modes_TE = [reshape((0:maxM).' + 0*(1:maxN), [], 1), ...
        reshape(0*(0:maxM).' + (1:maxN), [], 1)];
    modes_TM = [0, 0; modes_TE];
end

%% Get "waveguideMode" Objects
modesAll = [modes_TE; modes_TM];
modeTypes = [repmat("TE", size(modes_TE, 1), 1); ...
    repmat("TM", size(modes_TM, 1), 1)];

isRotated = false(size(modesAll, 1), 1);
modeTypes = [modeTypes; modeTypes(modesAll(:, 1) ~= 0)];
modesAll = [modesAll; modesAll(modesAll(:, 1) ~= 0, :)];
isRotated = [isRotated; true(size(modesAll, 1) - numel(isRotated), 1)];

%#ok<*AGROW>
for ii = 1:size(modesAll, 1)
    m = modesAll(ii, 1);
    n = modesAll(ii, 2);
    modeStructs(1, ii) = nLayer.getCoaxialModeStruct(...
        m, n, wgRi, wgRo, modeTypes(ii), isRotated(ii));
end

%% Sort by Cutoff
[~, sortInd] = sort([modeStructs.CutoffWavenumber]);
modeStructs = modeStructs(sortInd);

end

