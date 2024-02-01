clc;
clear;
close all;

%% Inputs
er = {4 - 0.4j};
ur = {1};
thk = {1.2};

f = linspace(32, 40, 2001);
% f = 36.61;

wgR = 30/128 * 25.4;
numModes = 3;

%% nLayerRectangularOld
tic;
NL1 = nLayerCircularTE_old(numModes, waveguideR=wgR, ...
    convergenceAbsTol=1e-6, verbosity=1);
toc;

%% nLayerRectangularFast
tic;
NL2 = nLayerCircularTE(numModes, waveguideR=wgR, ...
    verbosity=1);
toc;

%% Calculate
tic;
gam1 = NL1.calculate(f, er, ur, thk);
toc;

tic;
gam2 = NL2.calculate(f, er, ur, thk);
toc;

relErr = abs(max(gam1 - gam2))

%% Plot
figure;
plot(gam1, "-", LineWidth=1.5);
hold on;
plot(gam2, "-", LineWidth=1.5);
zplane([]);
grid on;
legend(["Original", "New"]);


% nLayer.plotModeStruct(NL2.modeStructs);

% figure;
% nLayerViewer(4, 1, 2, NL2, f);

