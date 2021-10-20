function [a, b] = setWaveguideBand(O, band, unitScaleFactor)
%SETWAVEGUIDEBAND Set waveguide a and b dimensions (mm) to a specific band
%   Sets the 

arguments
    O;
    band {mustBeText};
    unitScaleFactor(1, 1) {mustBeNumeric} = 1;
end

%% List of rectangular waveguide bands and dimensions in inches
bandNames = ["r",  "s",  "g",   "j",   "x", "ku",  "k",  "ka", "q",   "u",   "v",   "e",   "w" ];
bandDimsA = [4.30, 2.84, 1.872, 1.372, 0.9, 0.622, 0.42, 0.28, 0.224, 0.188, 0.148, 0.122, 0.10];
bandDimsB = [2.15, 1.34, 0.872, 0.622, 0.4, 0.311, 0.17, 0.14, 0.112, 0.094, 0.074, 0.061, 0.05];

%% Find specific band
bandIndex = find(strcmpi(bandNames, band));
if isempty(bandIndex)
    error("Rectangular waveguide band '%s' not found.", band);
end

%% Set dimensions
O.a = bandDimsA(bandIndex) .* 25.4 .* unitScaleFactor;
O.b = bandDimsB(bandIndex) .* 25.4 .* unitScaleFactor;

a = O.a;
b = O.b;

end

