clc; 
clear;
% close all;

%% Inputs
rPlot(:, 1) = linspace(-5, 40, 10001);
iPlot(1, :) = linspace(-20, 20, 10001);


x = 1;
n = 1;

Lc = 1;
L = 1;

%% Functions
alpha = @(x) x + 1j .* Lc .* x ./ (Lc + x);
alphaP = @(x) 1 + 1j .* Lc.^2 ./ (x + Lc).^2;

alpha = @(x) x + 1j .* Lc .* x ./ sqrt(Lc + x.^2);

alpha = @(x) x + 1j .* atan(x);

alpha = @(x) x + 1j .* Lc .* x ./ (Lc + abs(x));

alpha = @(x) x + 1j .* 2*atan(tanh(0.5 * abs(x)));

%% Calculate
val = alpha(rPlot + 1j*iPlot);

%% Plotting
figure;
showImage(rPlot, iPlot, -val, DisplayFormat="Phase");
clim([-180, 180]);
axis normal;
grid on;

% figure;
% showImage(rPlot, iPlot, val, DisplayFormat="dB");
% axis normal;
% grid on;





