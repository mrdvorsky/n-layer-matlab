clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(0, 2, ...
    modeSymmetryX="PEC", ...
    modeSymmetryY="PMC", ...
    modeSymmetryAxial="None", ...
    waveguideBand="Ka_TE01");

%% Plotting
for ii = 1:NL.numModes
    figure(Position=[100, 100, 1000, 800]);
    NL.waveguideModes(ii).showMode();
end






