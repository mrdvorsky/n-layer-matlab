function [modes_TE, modes_TM] = setModes(O, maxM, maxN)
%SETMODES Update modes in nLayerRectangular object and return a list of modes.
% Modes that don't satisfy symmetry conditions will not be included.
%
% Example Usage:
%   NL = nLayerRectangular(...);
%   [modesTE, modesTM] = NL.setModes(3, 2);     % 6 modes considered.
%
% Inputs:
%   maxM - maximum mode index m for any considered TEmn and TMmn modes.
%   maxN - maximum mode index n for any considered TEmn and TMmn modes.
% Outputs:
%   modes_TE - Rows of [m, n] mode index pairs for all considered TE modes.
%   modes_TM - Rows of [m, n] mode index pairs for all considered TM modes.
%
% Author: Matt Dvorsky

arguments
    O;
    maxM(1, 1) {mustBeInteger, mustBePositive} = 1;
    maxN(1, 1) {mustBeInteger, mustBeNonnegative} = 0;
end

%% Generate List of Modes
O.modes_TE = [reshape((1:2:maxM).' + 0*(0:2:maxN), [], 1), ...
    reshape(0*(1:2:maxM).' + (0:2:maxN), [], 1)];
O.modes_TM = O.modes_TE(find(O.modes_TE(:, 2) > 0), :);

%% Set Output
modes_TE = O.modes_TE;
modes_TM = O.modes_TM;

end


