function [modes_TE, modes_TM] = setModes(O, maxM, maxN, options)
%SETMODES Update modes in nLayerRectangular object and return a list of modes.
% Modes that don't satisfy symmetry conditions will not be included.
%
% Example Usage:
%   NL = nLayerRectangular(...);
%   [modesTE, modesTM] = NL.setModes(3, 2);     % 6 modes considered.
%   [modesTE, modesTM] = NL.setModes(3, 2, ModeSymmetryX="Odd");
%
% Inputs:
%   maxM - maximum mode index m for any considered TEmn and TMmn modes.
%   maxN - maximum mode index n for any considered TEmn and TMmn modes.
%
% Outputs:
%   modes_TE - Rows of [m, n] mode index pairs for all considered TE modes.
%   modes_TM - Rows of [m, n] mode index pairs for all considered TM modes.
%
% Named Arguments:
%   ModeSymmetryX ("even") - Symmetry condition about the x-axis.
%   ModeSymmetryY ("odd") - Symmetry condition about the y-axis.
%
% Author: Matt Dvorsky

arguments
    O;
    maxM(1, 1) {mustBeInteger, mustBeNonnegative} = 1;
    maxN(1, 1) {mustBeInteger, mustBeNonnegative} = 0;

    options.ModeSymmetryX {mustBeMember(options.ModeSymmetryX, ...
        ["None", "Even", "Odd"])} = O.modeSymmetryX;
    options.ModeSymmetryY {mustBeMember(options.ModeSymmetryY, ...
        ["None", "Even", "Odd"])} = O.modeSymmetryY;
end

%% Generate List of Modes
O.modes_TE = [reshape((0:maxM).' + 0*(0:maxN), [], 1), ...
    reshape(0*(0:maxM).' + (0:maxN), [], 1)];
O.modes_TE = O.modes_TE(2:end, :);  % Eliminate TE00 mode.

%% Eliminate Modes that Don't Satisfy Symmetry Conditions
if strcmp(options.ModeSymmetryX, "Even")
    O.modes_TE = O.modes_TE(mod(O.modes_TE(:, 1), 2) == 1, :);
elseif strcmp(options.ModeSymmetryX, "Odd")
    O.modes_TE = O.modes_TE(mod(O.modes_TE(:, 1), 2) == 0, :);
end

if strcmp(options.ModeSymmetryY, "Even")
    O.modes_TE = O.modes_TE(mod(O.modes_TE(:, 2), 2) == 1, :);
elseif strcmp(options.ModeSymmetryY, "Odd")
    O.modes_TE = O.modes_TE(mod(O.modes_TE(:, 2), 2) == 0, :);
end

%% Set TM modes.
O.modes_TM = O.modes_TE(O.modes_TE(:, 1) > 0 & O.modes_TE(:, 2) > 0, :);

%% Set Output
modes_TE = O.modes_TE;
modes_TM = O.modes_TM;

end


