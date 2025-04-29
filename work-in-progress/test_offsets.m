clc;
clear;
close all;

%% Inputs
er = {4 - 0.4j};
ur = {1};
thk = {1.2};

f = linspace(32, 40, 2001);

numModes = 3;

%% nLayer
NL1 = nLayerCircular(0, numModes, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
NL2 = nLayerCircular(0, numModes, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");

offsetX = num2cell(2 * ones(NL2.numModes, 1));
offsetY = num2cell(0 * ones(NL2.numModes, 1));
[NL2.modeStructs.OffsetX] = offsetX{:};
[NL2.modeStructs.OffsetY] = offsetY{:};

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



