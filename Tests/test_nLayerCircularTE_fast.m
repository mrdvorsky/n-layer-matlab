clc;
clear;
close all;

%% Inputs
er = {4 - 0.4j};
ur = {1};
thk = {1.2};

f = linspace(32, 40, 201);
% f = 36.61;

wgR = 30/128 * 25.4;
numModes = 5;

%% nLayerCircularTE_old
tic;
NL1 = nLayerCircularTE_old(numModes, waveguideR=wgR, ...
    convergenceAbsTol=1e-5, verbosity=1);
toc;

%% nLayerCircular
tic;
NL2 = nLayerCircular(0, numModes, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
NL2.calculate(f(1), er, ur, thk);
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

