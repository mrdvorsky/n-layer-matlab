clc;
clear;
close all;

%% Inputs
NL = nLayerRectangular(1, 0, ...
    modeSymmetryX="PEC", ...
    modeSymmetryY="PMC", ...
    modeSymmetryAxial="None", ...
    waveguideBand="Ka");

%% Plotting
for ii = flip(1:numel(NL.modeStructs))
    figure(Position=[100, 100, 1000, 800]);
    nLayer.plotModeStruct(NL.modeStructs(ii));
end


