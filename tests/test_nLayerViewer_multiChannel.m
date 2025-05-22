% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
er = {4 - 0.001j, 2 - 0.001j};
ur = {};
thk = {2, 0.1};

%% Create nLayer Object
NL = nLayerFilledRectangular(1, 0, waveguideBand="Ka");

%% Plot
figure;
nLayerViewer(er, ur, thk, NL, [26.5, 40]);

