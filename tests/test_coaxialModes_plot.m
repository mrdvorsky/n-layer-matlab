clc;
clear;
close all;

%% Inputs
NL = nLayerCoaxial(0, 2, ...
    modeSymmetryX="PEC", ...
    modeSymmetryY="PEC", ...
    modeSymmetryAxial="None", ...
    waveguideRi=0.3, ...
    waveguideRo=2.5);

%% Plotting
for ii = 1:NL.numModes
    figure(Position=[100, 100, 1000, 800]);
    NL.waveguideModes(ii).showMode();
end



