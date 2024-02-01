clc;
clear;
close all;

%% Inputs
er = {2.1 - 0.0001j};
ur = {1};
thk = {3};

f = linspace(26.5, 40, 801);

%% nLayerRectangularOld
tic;
NL1 = nLayerRectangular_old(3, 2, waveguideBand="ka", ...
    convergenceAbsTol=1e-5, verbosity=1);
toc;

%% nLayerRectangularFast
tic;
NL2 = nLayerRectangular(3, 2, waveguideBand="ka", ...
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


% nLayer.plotModeStruct(NL2.modeStructs)






