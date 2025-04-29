clc;
clear;
close all;

%% Inputs
f(:, 1) = linspace(8.2, 12.4, 800);

er = 2.1 - 0.01j;
ur = 1;
thk = 10;

NL = nLayerRectangular(3, 2, waveguideBand="X", modeSymmetryX="PEC", modeSymmetryY="PEC");
NL = nLayerRectangular(3, 2, waveguideBand="X");
% gam = NL.calculate(f, er, ur, thk);

%% Mode Cross Products
modeC = nLayer.modeCrossProduct(NL.waveguideModes(1), NL.waveguideModes(1));
modeC

%% Plot
% figure;
% plot(gam, "", LineWidth=1.5);
% hold on;
% zplane([]);


