clc;
clear;
close all;

%% Inputs
wgR = 5.8;
f = linspace(32, 40, 51).';

% Number of TE0n modes to consider
numModes = 3;

%% Create nLayer Object
NL = nLayerCircularTE(numModes, waveguideR=wgR, ...
    waveguideEr=1, waveguideUr=1);

NL2 = nLayerCircularTE(numModes, waveguideR=wgR/5, ...
    waveguideEr=25, waveguideUr=1);

%% Calculate structure 1
er = 25 - 0.01j;
thk = 0.7;

NL.printStructure(er, [], thk, Title="Case 1");
tic;
gam = NL.calculate(f, er, [], thk);
toc;

% Plot
figure;
plot(gam, ".-", Linewidth=1);
hold on;
title("Case 1");
zplane([]);
grid on;

