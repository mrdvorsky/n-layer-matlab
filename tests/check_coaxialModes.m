% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
NL = nLayerCoaxial(2, 2, waveguideRi=1.1, waveguideRo=2, ...
    modeSymmetryX="None", ...
    modeSymmetryY="None", ...
    modeSymmetryAxial="None");

%% Validate Symmetries
for ii = 1:NL.numModes
    NL.waveguideModes(ii).validateModeSymmetry();
    NL.waveguideModes(ii).validateModeAmplitude(...
        NumIntegralPointsPhi=8, ...
        NumIntegralPointsRho=16*16000, AbsTol=1e-4);
end

%% Plot
for ii = flip(1:NL.numModes)
    figure(Position=[300, 150, 1200, 800]);
    NL.waveguideModes(ii).showMode();
    title(sprintf("%s: (x='%s', y='%s', axial='%s')", ...
        NL.waveguideModes(ii).modeLabel, ...
        NL.waveguideModes(ii).symmetryX, ...
        NL.waveguideModes(ii).symmetryY, ...
        NL.waveguideModes(ii).symmetryAxial));
end

