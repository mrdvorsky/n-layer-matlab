clc;
clear;
close all;

%% Inputs
NL1 = nLayerCircular(0, 5, waveguideBand="Ka_TE01", modeSymmetryAxial="TE", ...
    integral_pointsKrc=250, ...
    integral_pointsKr=8192, ...
    integral_pointsPhi=64);
NL2 = nLayerCircular(0, 5, waveguideBand="Ka_TE01", modeSymmetryAxial="TE", ...
    integral_pointsKrc=200, ...
    integral_pointsKr=8192, ...
    integral_pointsPhi=64);

f = linspace(32, 40, 4001);
% f = 37.4751;
er = 4 - 0.0001j;
ur = 1;
thk = 2;

%% Calculate
gam1 = NL1.calculate(f, er, ur, thk);
gam2 = NL2.calculate(f, er, ur, thk);

%% Plot Error
figure;
plot(f, db(gam1 - gam2), "", LineWidth=1.5);
ylim([-300, 0]);
grid on;

%% PLot Values
figure;
plot(gam1, "", LineWidth=1.5);
hold on;
zplane([]);

