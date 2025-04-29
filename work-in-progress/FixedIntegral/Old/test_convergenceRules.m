clc;
clear;
close all;

%% Inputs
f = @(x) x ./ (1 + (0.1 - x).^2);
g = @(x) exp(-x.^2) .* exp(10j .* (x - 0.1).^2);

orderN = 120;
L = 10;

%% Integral
I1 = integral(@(x) f(x) .* g(x), -inf, inf)

[x, xw] = clenshawCurtisHalfOpen2(orderN, L, WeightingFunction=g);
I2 = sum(f(x) .* xw)

err_db = db(I2 - I1)

%% Plotting
xp = linspace(0, 10, 1001);

figure;
plot(xp, real(f(xp) .* g(xp)), "", LineWidth=1.5);
hold on;
plot(xp, imag(f(xp) .* g(xp)), "", LineWidth=1.5);
grid on;

figure;
plot(1:numel(x), f(x));


