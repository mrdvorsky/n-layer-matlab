function [modes] = getAllRectangularModes(m, n, wgA, wgB, options)
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
%   [modes] = getAllRectangularModes(m, n, a, b);
%
%   % Only return modes where the x-axis could be replaced with PEC.
%   [modes] = getAllRectangularModes(m, n, a, b, ...
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
        ["PEC", "PMC", "None"])} = "None";
    options.SymmetryY string {mustBeMember(options.SymmetryY, ...
        ["PEC", "PMC", "None"])} = "None";
    options.SymmetryAxial string {mustBeMember(options.SymmetryAxial, ...
        ["TE", "TM", "None"])} = "None";
end

%% Generate List of All Possible Modes
[TE_TM, m, n] = ndgrid(["TE", "TM"], unique(m), unique(n));

m = m(:);
n = n(:);
TE_TM = TE_TM(:);

%% Filter Out TE00, TMm0, and TM0n Modes
keepMode = (((m > 0) | (n > 0)) & TE_TM == "TE") ...
    | (((m > 0) & (n > 0)) & TE_TM == "TM");

%% Filter Modes by Symmetry
if strcmp(options.SymmetryY, "PMC")
    keepMode = keepMode & (mod(m, 2) == 1);
elseif strcmp(options.SymmetryY, "PEC")
    keepMode = keepMode & (mod(m, 2) == 0);
end

if strcmp(options.SymmetryX, "PMC")
    keepMode = keepMode & (mod(n, 2) == 1);
elseif strcmp(options.SymmetryX, "PEC")
    keepMode = keepMode & (mod(n, 2) == 0);
end

if ~strcmp(options.SymmetryAxial, "None")
    error("Axial mode symmetry not supported for rectangular waveguides.");
end

% Apply Filter
m = m(keepMode);
n = n(keepMode);
TE_TM = TE_TM(keepMode);

%% Get "waveguideMode" Objects
modes = waveguideMode.empty;
for ii = flip(1:numel(m))
    modes(1, ii) = waveguideMode.getRectangularMode(...
        m(ii), n(ii), wgA, wgB, TE_TM(ii));
end

end

