clc;
clear;
close all;

%% Inputs
er = {4 - 0.0j};
ur = {1};
thk = {10};

f = linspace(32, 40, 40001);
% f = 36.61;

wgR = 30/128 * 25.4;
numModes = 5;
f_res = 37.474;

%% nLayerCircularTE_old
tic;
NL1 = nLayerCircular(0, numModes, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
NL1.integral_pointsKrc = {200, 200, 200, 200, 200};
NL1.integral_pointsKr = {4096, 4*4096, 4*4096, 4*4096, 4*4096};
NL1.calculate(f_res, er, ur, thk);
toc;

%% nLayerCircular
tic;
NL2 = nLayerCircular(0, numModes, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
NL2.integral_pointsKrc = {50, 50, 150, 10, 10};
NL2.calculate(f_res, er, ur, thk);
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

figure;
plot(f, db(gam1 - gam2), "", LineWidth=1.5);
grid on;
ylim([-200, 0]);
xlim([min(f), max(f)]);


% nLayer.plotModeStruct(NL2.modeStructs);

% figure;
% nLayerViewer(4, 1, 2, NL2, f);

