clc;
clear;
close all;

%% Inputs
NL = nLayerCoaxial(2, 2, ...
    modeSymmetryX="None", ...
    modeSymmetryY="None", ...
    modeSymmetryAxial="None", ...
    waveguideRi=1.1, ...
    waveguideRo=2.5);

%% Plotting
for ii = flip(1:numel(NL.modeStructs))
    figure(Position=[100, 100, 1000, 800]);
    nLayer.plotModeStruct(NL.modeStructs(ii));
end


