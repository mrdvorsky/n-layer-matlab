clc;
clear;
close all;

%% Inputs
wgR = 15.8;

f = linspace(32, 40, 5001).';
numModes = 3;
convergenceTol = 1e-3;

%% Create nLayer Object
NL = nLayerCircularTE(numModes, waveguideR=wgR);

%% Calculate
er = [4 - 0.001j];
thk = [1.7];

Smn = NL.calculate(f, er, [], thk);

%% Plot
figure;
plot(Smn(:, 1), LineWidth=1.5);
hold on;
zplane([]);

