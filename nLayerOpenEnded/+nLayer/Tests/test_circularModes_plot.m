clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(0, 5, ...
    modeSymmetryAxial="TE", ...
    waveguideR=2.5);

%% Plotting
for ii = flip(1:numel(NL.waveguideModes))
    figure(Position=[100, 100, 1000, 800]);
    nLayer.plotModeStruct(NL.modeStructs(ii));
end
