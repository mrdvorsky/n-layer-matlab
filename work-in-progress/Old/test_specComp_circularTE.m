clc;
clear;
close all;

%% Inputs
wgR = 30/128*25.4;

[specEx, specEy, kc] = nLayer.getSpectrumCircular(wgR, 0, 1, "TE");

%% Get Integrand
modeFun1 = @(x) x .* besselj(1, wgR*x).^2 ...
    ./ (x.^2 - kc.^2).^2;

modeFun2 = @(x) x .* specEy(x, 0, x, 0).^2;

%% Plotting
krPlot = linspace(0.001, 5, 10000);

figure;
plot(krPlot, real(modeFun1(krPlot)), "", LineWidth=1.5);
grid on;

figure;
plot(krPlot, real(modeFun2(krPlot)), "", LineWidth=1.5);
grid on;




