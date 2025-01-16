clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(2, 2, ...
    modeSymmetryX="None", ...
    modeSymmetryY="None", ...
    modeSymmetryAxial="None", ...
    waveguideR=2.5);

%% Plotting
for ii = flip(1:numel(NL.modeStructs))
    figure(Position=[100, 100, 1000, 800]);
    nLayer.plotModeStruct(NL.modeStructs(ii));
end
