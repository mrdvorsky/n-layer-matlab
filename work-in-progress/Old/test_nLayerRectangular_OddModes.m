clc;
clear;
close all;

%% Inputs
er = [2.1 - 0.001j];
ur = [1];
thk = [10];

f = linspace(8.2, 12.4, 801);

%% nLayerRectangular
tic;
NL2 = nLayerRectangular(3, 2, waveguideBand="x", modeSymmetryX="None", ...
    convergenceAbsTol=1e-4, verbosity=1);
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

% figure;
% plot(f, real(gam2), "-", LineWidth=1.5);
% hold on;
% plot(f, imag(gam2), "-", LineWidth=1.5);
% grid on;
% legend(["Real", "Imag"]);






