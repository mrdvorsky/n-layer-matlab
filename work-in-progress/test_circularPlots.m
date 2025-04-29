clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(2, 2, waveguideBand="Ka_TE01", modeSymmetryX="None", modeSymmetryY="None");
% NL = nLayerRectangular(2, 2, waveguideBand="W", modeSymmetryX="None", modeSymmetryY="None");

% NL.calculate(100, {1}, {2}, {2});

cellTest = num2cell(ones(NL.numModes, 1));
[NL.modeStructs.OffsetX] = cellTest{:};

%% Plotting
nLayer.plotModeStruct(NL.modeStructs);


