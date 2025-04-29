clc;
clear;
close all;

%% Inputs
f = linspace(32, 40, 2001);

er = [4 - 0.0004j, 2.1 - 0.001j];
thk = [1, 0.25];

%% nLayer Object
NL = nLayerCircularTE(3, waveguideR=29/128 * 25.4);

%% Calculate
% gam = NL.calculate(f, er, [], thk);

S11plusS21 = NL.calculate(f, [er, 1], [1, 1, 1 - 1j*1e15], [thk, inf]);
S11minusS21 = NL.calculate(f, er, [], thk);

S11 = 0.5 * (S11plusS21 + S11minusS21);
S21 = 0.5 * (S11plusS21 - S11minusS21);

%% Plot
% figure;
% p = plot(gam, "-", LineWidth=1.5);
% hold on;
% zplane([]);

figure;
plot(S11, "-", LineWidth=1.5);
hold on;
plot(S21, "-", LineWidth=1.5);
zplane([]);

figure;
plot(S11plusS21, "-", LineWidth=1.5);
hold on;
plot(S11minusS21, "-", LineWidth=1.5);
zplane([]);



