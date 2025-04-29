clc;
clear;
close all;

%% Inputs
er = {2.1 - 0.0005j};
ur = {1};
thk = {15};

f = linspace(8.2, 12.4, 801);

%% nLayerRectangularOld
tic;
NL1 = nLayerRectangular_old(3, 2, waveguideBand="X", ...
    convergenceAbsTol=1e-6, verbosity=1);
toc;

%% nLayerRectangularFast
tic;
NL2 = nLayerRectangular(3, 2, waveguideBand="X");
NL2.calculate(f(1), er, ur, thk);
toc;

%% Calculate
tic;
gam1 = NL1.calculate(f, er, ur, thk);
toc;

tic;
for ii = 1:1
    gam2 = NL2.calculate(f, er, ur, thk);
end
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


% NL3 = nLayerRectangular(1, 0, waveguideBand="ka");

% figure;
% nLayerViewer([1, 4], [], [1, 1], NL2, f, NL3, f, NumFrequencySamplesPerMarker=10);

% nLayer.plotModeStruct(NL2.modeStructs)






