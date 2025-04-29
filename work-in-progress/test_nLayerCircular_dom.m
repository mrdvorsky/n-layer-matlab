clc;
clear;
close all;

%% Inputs
f = linspace(33, 40, 800);

er = {2.1 - 0.002j};
thk = {3};
ur = {1};

%% nLayer Object
tic;
NL = nLayerCircular(1, 3, waveguideBand="Ka_High");
NL.calculate(f(1), er, ur, thk);
toc;

%% Calculate
NL.receiveModeIndices = 1:NL.numModes;
gam = NL.calculate(f, er, ur, thk);


figure;
plot(gam, "", LineWidth=1.5);
hold on;
zplane([]);

figure;
plot(f, db(gam), "", LineWidth=1.5);
grid on;
legend(NL.modeLabels);

% figure;
% nLayerViewer(er, ur, thk, NL, f);





