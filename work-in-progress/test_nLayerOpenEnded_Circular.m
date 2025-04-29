clc;
clear;
close all;

%% Inputs
er = [2.1 - 0.01j];
ur = [1];
thk = [10];

f = linspace(10, 12.4, 801);

%% nLayerCircular
tic;
NL1 = nLayerCircular(3, 2, waveguideR=30, ...
    convergenceAbsTol=1e-4, verbosity=1);
toc;

%% Calculate
tic;
gam1 = NL1.calculate(f, er, ur, thk);
toc;

%% Plot
figure;
plot(gam1, "-", LineWidth=1.5);
hold on;
% plot(gam2, "-", LineWidth=1.5);
zplane([]);
grid on;
% legend(["Original", "New"]);

% figure;
% plot(f, real(gam2), "-", LineWidth=1.5);
% hold on;
% plot(f, imag(gam2), "-", LineWidth=1.5);
% grid on;
% legend(["Real", "Imag"]);






