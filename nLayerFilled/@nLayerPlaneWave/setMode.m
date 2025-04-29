function [] = setMode(O, modeTE_m, modeTE_n)
%SETMODES Update mode in a nLayerFilledRectangular object.
%
% Example Usage:
%   NL = nLayerFilledRectangular(...);
%   [modesTE, modesTM] = NL.setMode(2, 1);      % Use TE21 mode.
%
% Inputs:
%   modeTE_m - Mode index m for the considered TEmn mode.
%   modeTE_n - Mode index n for the considered TEmn mode.
%
% Author: Trent Moritz

arguments
    O;
    modeTE_m(1, 1) {mustBeInteger, mustBeNonnegative} = 1;
    modeTE_n(1, 1) {mustBeInteger, mustBeNonnegative} = 0;
end

%% Set Mode
O.modeTE_m = modeTE_m;
O.modeTE_n = modeTE_n;

end


