clc;
clear;
close all;

%% Inputs
rPlot(:, 1) = linspace(0, 40, 1001);
iPlot(1, :) = linspace(-20, 20, 1001);


x = 1;
n = 1;

Lc = 1;
L = 1;

%% Functions
alpha = @(x) x + 1j .* Lc .* x ./ (Lc + x);
alphaP = @(x) 1 + 1j .* Lc.^2 ./ (x + Lc).^2;

fun = @(k) exp(-1j*x .* k) .* (exp(-1j * x .* (alpha(k) - k)) .* alphaP(k) ...
    .* (cos(2*n * acot(sqrt(k ./ L)))) - exp(1));

%% Calculate
val = fun(rPlot + 1j*iPlot);

%% Plotting
figure;
showImage(rPlot, iPlot, -val, DisplayFormat="Phase");
axis normal;
grid on;

figure;
showImage(rPlot, iPlot, -val, DisplayFormat="dB");
axis normal;
grid on;





