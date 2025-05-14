% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
f(:, 1) = linspace(32, 40, 800);

er = 4 - 0.001j;
ur = 1;
thk = 2;

%% Create Object
NL = nLayerCircular(0, 5, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");

%% Calculate
nLayer.printStructure(er, ur, thk);
gam = NL.calculate(f, er, ur, thk);

%% Plot
figure;
plotComplex(f, gam, ".-", LineWidth=1.5);
legend(NL.getOutputLabels());
grid on;

