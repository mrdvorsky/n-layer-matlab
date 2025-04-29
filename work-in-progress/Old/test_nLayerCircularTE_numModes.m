clc;
clear;
close all;

%% Inputs
wgR = 29/128 * 25.4;

f(:, 1) = linspace(32, 40, 2001);

% Number of TE0n modes to consider
numModes = 30;

% Tolerance for convergence (-60 dB)
convergenceTol = 1e-4;

%% Create nLayer Object
NL = nLayerCircularTE(numModes, waveguideR=wgR);
NL.convergenceAbsTol = convergenceTol;
NL.verbosity = 1;

%% Calculate Structure
er2 = 4 - 0.0001j;
ur2 = [];
thk2 = 2;

NL.printStructure(er2, ur2, thk2(1), Title="Case 2");
tic;
gam2 = NL.calculate(f, er2, ur2, thk2);
toc;

%% Plotting
figure;
plot(gam2, ".-", Linewidth=1);
hold on;
title("Case 2");
zplane([]);
grid on;

figure;
plot(f, rad2deg(angle(gam2)), "", LineWidth=1.5);
xlabel("f (GHz)");
ylabel("S_{11}");
grid on;






