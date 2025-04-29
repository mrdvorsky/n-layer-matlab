clc;
clear;
close all;

%% Inputs
wgR = 1;

xView = 1.1 * linspace(-1, 1, 1001);
yView = 1.1 * linspace(-1, 1, 1001);

m = 1;
n = 1;

%% Calculate Spectrums
[specEx{1}, specEy{1}, modeCutoffs, phaseScale] ...
    = nLayerOpenEnded.getSpectrumCircular(wgR, m, n, "TE");

modeStruct = nLayerOpenEnded.createModeStruct(...
    SpecEx_TE=specEx, ...
    SpecEy_TE=specEy, ...
    CutoffBeta_TE=modeCutoffs, ...
    PhaseScaleFactor_TE=phaseScale, ...
    OutputModes_TE=true);

nLayerOpenEnded.plotModeStruct(modeStruct, ...
    SizeX=2, SizeY=2, PlotCoordinates="Rectangular");




