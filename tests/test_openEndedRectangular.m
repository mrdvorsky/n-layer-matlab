% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
f(:, 1) = linspace(8.2, 12.4, 800);

er = 2.1 - 0.001j;
ur = 1;
thk = 10;

%% Create Object
NL = nLayerRectangular(3, 2, waveguideBand="X");

%% Calculate
nLayer.printStructure(er, ur, thk);
gam = NL.calculate(f, er, ur, thk);

%% Plot
figure;
plotComplex(f, gam, "", LineWidth=1.5);
legend(NL.getOutputLabels());
grid on;

