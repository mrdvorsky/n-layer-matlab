clc;
clear;
close all;

%% Inputs
f(:, 1) = linspace(32, 40, 800);

er = 4 - 0.001j;
ur = 1;
thk = 2;

NL = nLayerCircular(0, 5, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
gam = NL.calculate(f, er, ur, thk);

%% Plot
figure;
plot(gam, ".-", LineWidth=1.5);
hold on;
zplane([]);


