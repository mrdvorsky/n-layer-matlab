clc;
clear;
close all;

%% Inputs
er = [2.1 - 0.01j];
ur = [1];
thk = [10];

f = linspace(8.2, 12.4, 801);

%% nLayerRectangularOld
% tic;
% NL1 = nLayerRectangularOld(3, 2, waveguideBand="x", ...
%     convergenceAbsTol=1e-4, verbosity=1);
% toc;

%% nLayerRectangular
tic;
NL2 = nLayerRectangular(3, 2, waveguideBand="x", ...
    convergenceAbsTol=1e-4, verbosity=1, ...
    modeSymmetryY="Odd", modeSymmetryX="Even");
toc;
NL2.waveguideEr = 1;

%% Calculate
% tic;
% gam1 = NL1.calculate(f, er, ur, thk);
% toc;

tic;
gam2 = NL2.calculate(f, er, ur, thk);
toc;

relErr = abs(max(gam1 - gam2))

%% Plot
figure;
% plot(gam1, "-", LineWidth=1.5);
hold on;
plot(gam2, "-", LineWidth=1.5);
zplane([]);
grid on;
legend(["Original", "New"]);

figure;
plot(f, real(gam2), "-", LineWidth=1.5);
hold on;
plot(f, imag(gam2), "-", LineWidth=1.5);
grid on;
legend(["Real", "Imag"]);






