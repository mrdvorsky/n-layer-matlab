clc;
clear;
close all;

%% Inputs
N = 32000;
Nk = 16000;

k = 800 * linspace(-1, 1, Nk);
xPlot = linspace(-2.2, 2.2, 10001);

%% Calculate Spectrum
[x, w] = fejer2(N, -1, 1);

spec1 = nufft(w, x ./ 2*pi, k);

%% Calculate Conv
convX = nufft(spec1.^1, -k ./ 2*pi, xPlot) .* abs(k(2) - k(1));

%% Plotting
figure;
plot(xPlot, real(convX), "", LineWidth=1.5);
grid on;

% figure;
% plot(xPlot, real(convX) - , "", LineWidth=1.5);
% grid on;






