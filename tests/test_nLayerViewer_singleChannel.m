clc;
clear;
close all;

%% Inputs
er = {4 - 0.001j, 2 - 0.001j};
ur = {};
thk = {2, 0.1};

NL = nLayerCircular(0, 5, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
% NL = nLayerRectangular(5, 4, waveguideBand="Ka");

%% Plot
figure;
nLayerViewer(er, ur, thk, NL, [32, 40]);


