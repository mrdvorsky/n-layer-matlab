clc;
clear;
close all;

%% Inputs
NL = nLayerCoaxial(2, 2, waveguideRi=0.2, waveguideRo=2, ...
    modeSymmetryX="None", ...
    modeSymmetryY="None", ...
    modeSymmetryAxial="None");

%% Validate Symmetries
for ii = 1:NL.numModes
    NL.waveguideModes(ii).validateModeSymmetry();
end

%% Plot
for ii = flip(1:NL.numModes)
    figure;
    NL.waveguideModes(ii).showMode();
    title(sprintf("%s: (x='%s', y='%s', axial='%s')", ...
        NL.waveguideModes(ii).modeLabel, ...
        NL.waveguideModes(ii).symmetryX, ...
        NL.waveguideModes(ii).symmetryY, ...
        NL.waveguideModes(ii).symmetryAxial));
end

