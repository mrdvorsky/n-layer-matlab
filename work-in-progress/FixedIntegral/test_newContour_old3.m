clc;
clear;
close all;

%% Inputs
f = @(x) 1 ./ (1 + x).^2;

Lc = 1;

%% Helper
krToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + (x));
krToKrcPrime = @(x) 1 + 1j .* Lc.^2 ./ ((x) + Lc).^2;

krcToKr1 = @(x) 0.5 * (sqrt(x.^2 + Lc*(2-2j)*x + 2j*Lc.^2) + x - Lc - 1j*Lc);
krcToKr2 = @(x) 0.5 * (-sqrt(x.^2 + Lc*(2-2j)*x + 2j*Lc.^2) + x - Lc - 1j*Lc);

switchKrc = @(x) angle(1j * (x.^2 + 2*x + 1 + 1j) ./ (1 + x)) > 0;

krcToKr = @(x) krcToKr1(x) .* (switchKrc(x)) + krcToKr2(x) .* (~switchKrc(x));

%% Integral
% I1 = integral(f, 0, inf);
% I2 = integral(@(x) f(krToKrc(x)) .* krToKrcPrime(x), 0, inf);
% I3 = integral(@(x) f(krcToKr(x)) .* krcToKrPrime(x), 0, inf);

%% Plotting
rPlot(:, 1) = linspace(-10, 10, 1001);
iPlot(1, :) = linspace(-10, 10, 1001);
cPlot = rPlot + 1j*iPlot;

figure;
showImage(rPlot, iPlot, krToKrc(cPlot), DisplayFormat="Phase");
axis normal;
colormap hsv;

figure;
showImage(rPlot, iPlot, krToKrc(cPlot), DisplayFormat="Magnitude");
axis normal;

% figure;
% showImage(rPlot, iPlot, cPlot - krcToKr(krToKrc(cPlot)), DisplayFormat="dB");
% axis normal;






