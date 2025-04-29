clc;
clear;
close all;

%% Inputs
thetaR(:, 1) = linspace(0, pi/2, 101);

%% Calculate
thetaR = thetaR(2:end-1);
xR = cot(thetaR);

x = xR + 1j*xR ./ (1 + abs(xR));
theta = acot(x);

theta = thetaR + 1j .* thetaR.^2 .* cos(thetaR).^2;
x = cot(theta);

%% Plotting
figure;
plot(thetaR + 0.000001j, "", LineWidth=1.5);
grid on;
hold on;
zplane([]);

figure;
plot(theta, "", LineWidth=1.5);
grid on;
hold on;
zplane([]);

figure;
plot(x, "", LineWidth=1.5);
grid on;
hold on;
zplane([]);