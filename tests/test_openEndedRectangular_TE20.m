% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
f(:, 1) = linspace(32, 40, 800);

er = 4 - 0.1j;
ur = 1;
thk = 2;

%% Create Object
NL = nLayerRectangular(5, 4, waveguideBand="Ka", modeSymmetryX="PEC", modeSymmetryY="PEC");
NL.waveguideA = 1.4 * NL.waveguideA;

NL.frequencyRange = f;

%% Calculate
nLayer.printStructure(er, ur, thk);
gam = NL.calculate(f, er, ur, thk);

%% Plot
figure;
plotComplex(f, gam, ".-", LineWidth=1.5);
legend(NL.getOutputLabels());
grid on;
