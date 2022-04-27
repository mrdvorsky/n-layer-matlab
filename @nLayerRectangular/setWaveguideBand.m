function [a, b] = setWaveguideBand(O, band, options)
%SETWAVEGUIDEBAND Set waveguide a and b dimensions (mm) to a specific band.
% Calling this functions sets the broad and narrow dimensions, O.a and O.b,
% to the dimensions of the band specified. The default unit is mm, however,
% this can be changed by specifying a the UnitScaleFactor optional
% parameter, which will be used as a multiplier to determine the new values
% of O.a and O.b.
%
% After calling this function, the "recomputeInterpolants" function should
% be called before calling "calculate".
%
% Example Usage:
%   NL = nLayerRectangular(...);
%   NL.setWaveguideBand(band, UnitScaleFactor=0.001);
%   NL.recomputeInterpolants();
%
% Inputs:
%   band - Case-Insensitive waveguide band designator.
% Outputs:
%   a - New value of O.a (waveguide broad dimension).
%   b - New value of O.b (waveguide narrow dimension).
% Named Arguments:
%   UnitScaleFactor (1) - Multiplier for O.a and O.b.
%
% Author: Matt Dvorsky

arguments
    O;
    band {mustBeTextScalar};
    options.UnitScaleFactor(1, 1) {mustBeReal, mustBePositive} = 1;
end

%% List of rectangular waveguide bands and dimensions in inches
bandNames = ["r",  "s",  "g",   "j",   "x", "ku",  "k",  "ka", "q",   "u",   "v",   "e",   "w" ];
bandDimsA = [4.30, 2.84, 1.872, 1.372, 0.9, 0.622, 0.42, 0.28, 0.224, 0.188, 0.148, 0.122, 0.10];
bandDimsB = [2.15, 1.34, 0.872, 0.622, 0.4, 0.311, 0.17, 0.14, 0.112, 0.094, 0.074, 0.061, 0.05];

%% Find specific band
bandIndex = find(strcmpi(bandNames, band));
if isempty(bandIndex)
    error(strcat("Rectangular waveguide band '%s' not found. Must ", ...
        "match one of the following (case insensitive): {'%s'}."), band, ...
        strjoin(bandNames, "', '"));
end

%% Set dimensions
O.a = bandDimsA(bandIndex) .* 25.4 .* options.UnitScaleFactor;
O.b = bandDimsB(bandIndex) .* 25.4 .* options.UnitScaleFactor;

a = O.a;
b = O.b;

end

