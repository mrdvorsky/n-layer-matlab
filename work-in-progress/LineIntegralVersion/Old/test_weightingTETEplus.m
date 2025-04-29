clc;
clear;
close all;

%% Inputs
Nkr = 10000;

krMax = 500;

r(:, 1) = 10.^linspace(-4, 1, 1000);

%% Integral
% [kr(1, :), kr_w(1, :)] = fejer2_halfOpen(Nkr, 1);
[kr(1, :), kr_w(1, :)] = fejer2(Nkr, 0, krMax);

vals = sum(besselj(2, kr .* r) .* kr_w, 2);

%% Plotting
figure;
semilogx(r, abs(vals), "", LineWidth=1.5);
grid on;












