function [modes] = getAllCoaxialModes(m, n, wgRi, wgRo, options)
%Get all possible "waveguideMode" objects for a coaxial waveguide.
% This functions returns "waveguideMode" objects for all coaxial
% waveguide modes that match the pattern TEmn or TMmn, for all
% combinations of the input vectors "m" and "n";
%
% Optionally, symmetry filters can be applied, so that all returned mode
% objects have the specified symmetry.
%
% Example Usage:
%   % All modes, regardless of symmetry.
%   [modes] = getAllCoaxialModes(m, n, wgRi, wgRo);
%
%   % Only return modes where the x-axis could be replaced with PEC.
%   [modes] = getAllCoaxialModes(m, n, wgRi, wgRo, ...
%               modeSymmetryX="PEC");
%
%
% Inputs:
%   m - Vector of "m" values for returned TEmn and TMmn modes.
%   n - Vector of "n" values for returned TEmn and TMmn modes.
%   wgRi - Inner radius of coaxial waveguide.
%   wgRo - Inner radius of coaxial waveguide.
%
% Author: Matt Dvorsky


arguments
    m(:, 1) {mustBeInteger, mustBeNonnegative};
    n(1, :) {mustBeInteger, mustBeNonnegative};
    wgRi(1, 1) {mustBePositive};
    wgRo(1, 1) {mustBePositive};

    options.SymmetryX string {mustBeMember(options.SymmetryX, ...
        ["PEC", "PMC", "None"])} = "None";
    options.SymmetryY string {mustBeMember(options.SymmetryY, ...
        ["PEC", "PMC", "None"])} = "None";
    options.SymmetryAxial string {mustBeMember(options.SymmetryAxial, ...
        ["TE", "TM", "None"])} = "None";
end

%% Generate List of All Possible Modes
[isRotated, TE_TM, n, m] = ndgrid([false, true], ["TE", "TM"], ...
    unique(n), unique(m));

m = m(:);
n = n(:);
TE_TM = TE_TM(:);
isRotated = isRotated(:);

%% Filter Out Rotated TE0n and TM0n Modes
keepMode = ~((m == 0) & isRotated) & (n > 0) ...
    | ((m == 0) & (n == 0) & strcmp(TE_TM, "TM") & ~isRotated);

%% Filter Modes by Symmetry
if strcmp(options.SymmetryX, "PEC")
    keepMode = keepMode ...
        & xor(isRotated, strcmp(TE_TM, "TE"));
elseif strcmp(options.SymmetryX, "PMC")
    keepMode = keepMode ...
        & xor(~isRotated, strcmp(TE_TM, "TE"));
end

if strcmp(options.SymmetryY, "PEC")
    keepMode = keepMode ...
        & xor(xor((mod(m, 2) == 1), isRotated), strcmp(TE_TM, "TE"));
elseif strcmp(options.SymmetryY, "PMC")
    keepMode = keepMode ...
        & xor(xor((mod(m, 2) == 0), isRotated), strcmp(TE_TM, "TE"));
end

if strcmp(options.SymmetryAxial, "TE")
    keepMode = keepMode ...
        & ((m == 0) & strcmp(TE_TM, "TE"));
elseif strcmp(options.SymmetryAxial, "TM")
    keepMode = keepMode ...
        & ((m == 0) & strcmp(TE_TM, "TM"));
end

% Apply Filter
m = m(keepMode);
n = n(keepMode);
TE_TM = TE_TM(keepMode);
isRotated = isRotated(keepMode);

%% Get "waveguideMode" Objects
modes = waveguideMode.empty;
for ii = flip(1:numel(m))
    modes(1, ii) = waveguideMode.getCoaxialMode(...
        m(ii), n(ii), wgRi, wgRo, TE_TM(ii), isRotated(ii));
end

end

