function [modeStructs] = getRectangularModes(maxM, maxN, wgA, wgB, options)
%GETRECTANGULARMODES Get "waveguideMode" objects for rectangular waveguide.
% 
% Author: Matt Dvorsky

arguments
    maxM(1, 1) {mustBeInteger, mustBeNonnegative};
    maxN(1, 1) {mustBeInteger, mustBeNonnegative};
    wgA(1, 1) {mustBePositive};
    wgB(1, 1) {mustBePositive};

    options.SymmetryX string {mustBeMember(options.SymmetryX, ...
        ["PEC", "PMC", "None"])} = "PEC";
    options.SymmetryY string {mustBeMember(options.SymmetryY, ...
        ["PEC", "PMC", "None"])} = "PMC";
    options.SymmetryAxial string {mustBeMember(options.SymmetryAxial, ...
        ["TE", "TM", "None"])} = "None";
end

%% Generate List of Modes
modes_TE = [reshape((0:maxM).' + 0*(0:maxN), [], 1), ...
    reshape(0*(0:maxM).' + (0:maxN), [], 1)];
modes_TE = modes_TE(2:end, :);  % Eliminate TE00 mode.

%% Eliminate Modes that Don't Satisfy Symmetry Conditions
if strcmp(options.SymmetryY, "PMC")
    modes_TE = modes_TE(mod(modes_TE(:, 1), 2) == 1, :);
elseif strcmp(options.SymmetryY, "PEC")
    modes_TE = modes_TE(mod(modes_TE(:, 1), 2) == 0, :);
end

if strcmp(options.SymmetryX, "PMC")
    modes_TE = modes_TE(mod(modes_TE(:, 2), 2) == 1, :);
elseif strcmp(options.SymmetryX, "PEC")
    modes_TE = modes_TE(mod(modes_TE(:, 2), 2) == 0, :);
end

if ~strcmp(options.SymmetryAxial, "None")
    error("Axial mode symmetry not supported for rectangular waveguides.");
end

%% Set TM Modes
modes_TM = modes_TE(modes_TE(:, 1) > 0 & modes_TE(:, 2) > 0, :);

%% Get "waveguideMode" Objects
modesAll = [modes_TE; modes_TM];
modeTypes = [repmat("TE", size(modes_TE, 1), 1); ...
    repmat("TM", size(modes_TM, 1), 1)];

%#ok<*AGROW>
for ii = 1:size(modesAll, 1)
    m = modesAll(ii, 1);
    n = modesAll(ii, 2);
    modeStructs(1, ii) = nLayer.getRectangularModeStruct(...
        m, n, wgA, wgB, modeTypes(ii));
end

%% Sort by Cutoff
[~, sortInd] = sort([modeStructs.CutoffWavenumber]);
modeStructs = modeStructs(sortInd);

end

