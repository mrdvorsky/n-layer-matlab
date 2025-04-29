clc;
clear;
close all;

%% Inputs
f(:, 1) = linspace(8.2, 12.4, 800);

er = 2.1 - 0.01j;
ur = 1;
thk = 10;

NL = nLayerFilledRectangular(1, 0, waveguideBand="X");

%% Plot
gam = NL.calculate(f, er, ur, thk);

figure;
plot(gam(:, :), "", LineWidth=1.5);
hold on;
zplane([]);


