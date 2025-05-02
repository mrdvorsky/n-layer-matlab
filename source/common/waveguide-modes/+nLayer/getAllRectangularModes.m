function [waveguideModes] = getAllRectangularModes(m, n, wgA, wgB, options)
%Get all possible "waveguideMode" objects for a rectangular waveguide.
% This functions returns "waveguideMode" objects for all rectangular
% waveguide modes that match the pattern TEmn or TMmn, for all
% combinations of the input vectors "m" and "n";
%
% Optionally, symmetry filters can be applied, so that all returned mode
% objects have the specified symmetry.
%
% Example Usage:
%   % All modes, regardless of symmetry.
%   [waveguideModes] = getAllRectangularModes(m, n, a, b);
%
%   % Only return modes where the x-axis could be replaced with PEC.
%   [waveguideModes] = getAllRectangularModes(m, n, a, b, ...
%               modeSymmetryX="PEC");
%
%
% Inputs:
%   m - Vector of "m" values for returned TEmn and TMmn modes.
%   n - Vector of "n" values for returned TEmn and TMmn modes.
%   wgA - Waveguide length along x-dimension.
%   wgB - Waveguide length along y-dimension.
%
% Author: Matt Dvorsky

arguments
    m(:, 1) {mustBeInteger, mustBeNonnegative};
    n(1, :) {mustBeInteger, mustBeNonnegative};
    wgA(1, 1) {mustBePositive};
    wgB(1, 1) {mustBePositive};

    options.SymmetryX string {mustBeMember(options.SymmetryX, ...
        ["PEC", "PMC", "None"])} = "PEC";
    options.SymmetryY string {mustBeMember(options.SymmetryY, ...
        ["PEC", "PMC", "None"])} = "PMC";
    options.SymmetryAxial string {mustBeMember(options.SymmetryAxial, ...
        ["TE", "TM", "None"])} = "None";
end

%% Generate List of All Possible Modes
rotated(1, 1, :) = [false, true];
[m, n, rotated] = broadcastArrays(unique(m), unique(n), rotated);
m = m(:);
n = n(:);
rotated = rotated(:);

%% Filter Modes by Symmetry
keepMode = true(size(m));
if strcmp(options.SymmetryY, "PMC")
    
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
waveguideModes = nLayer.waveguideMode.empty;
for ii = flip(1:size(modesAll, 1))
    m = modesAll(ii, 1);
    n = modesAll(ii, 2);
    waveguideModes(1, ii) = nLayer.getRectangularModeStruct(...
        m, n, wgA, wgB, modeTypes(ii));
end

%% Sort by Cutoff
[~, sortInd] = sort([waveguideModes.kc0]);
waveguideModes = waveguideModes(sortInd);

end

