clc;
clear;
close all;

%% Inputs
f(:, 1) = linspace(32, 40, 800);

er = 4 - 0.001j;
ur = 1;
thk = 2;

NL = nLayerCircular(0, 5, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");

NLC = nLayerClamped(NL);

%% Mode Cross Products
modeC = nLayer.modeCrossProduct(NL.waveguideModes(2), NL.waveguideModes(2));
modeC

NLsolver = nLayerInverse(2);

%% Plot
gam = NLC.calculate(f, er, ur, thk);

figure;
plot(gam(:, :), ".-", LineWidth=1.5);
hold on;
zplane([]);


