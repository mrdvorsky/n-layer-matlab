function [modesTE, modesTM] = setModes(O, maxM, maxN)
%SETMODES Update modes in nLayerRectangular object and return a list of modes.
% Inputs:
%   maxM - maximum mode index m for any considered TEmn and TMmn modes.
%   maxN - maximum mode index n for any considered TEmn and TMmn modes.
% Outputs:
%   modesTE - Rows of [m, n] mode index pairs for all considered TE modes.
%   modesTM - Rows of [m, n] mode index pairs for all considered TM modes.
%
% After calling this function, the "recomputeInterpolants" function should
% be called before calling "calculate".
%
% Example Usage:
%   NL = nLayerRectangular(...);
%   [modesTE, modesTM] = NL.setModes(3, 2);     % 6 modes considered.
%   NL.recomputeInterpolants();
%
% Author: Matt Dvorsky

arguments
    O;
    maxM(1, 1) = 1;
    maxN(1, 1) = 0;
end

%% Generate List of Modes
O.modesTE = [reshape((1:2:maxM).' + 0*(0:2:maxN), [], 1), ...
    reshape(0*(1:2:maxM).' + (0:2:maxN), [], 1)];
O.modesTM = O.modesTE(find(O.modesTE(:, 2) > 0), :);

O.numModes = size(O.modesTE, 1) + size(O.modesTM, 1);

%% Set Output
modesTE = O.modesTE;
modesTM = O.modesTM;

end


