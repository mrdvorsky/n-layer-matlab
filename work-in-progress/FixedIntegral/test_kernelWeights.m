clc;
clear;
close all;

%% Inputs
N = 1000;
L = 1;
Lc = 1;

k = 10;
n = 1;

%% Get Kernels
[krc, momentKernel] = getContourIntegrals(N, L, Lc);

w = @(x) momentKernel(x, k);
w2 = @(x) w(cot(0.5*x).^2) .* (1 + 0*sin(n*x));
w2 = @(x) w(cot(0.5*x).^2) .* (0 + 1*exp(-1j * n*x));
% w2 = @(x) w(x) .* (1 + 0*sin(n*x));
M_k = integral(w2, 0, pi)

%% Plotting
rPlot(:, 1) = linspace(-1, 4, 1000);
iPlot(1, :) = linspace(-2, 2, 1000);

figure;
showImage(rPlot, iPlot, w2(rPlot + 1j*iPlot), DisplayFormat="dB");
axis normal;










