clc;
clear;
close all;

%% Inputs
f = @(x) 1 ./ (1 + x).^2;


%% Helper
krToKrc = @(z) z .* exp(0.25j * pi ./ (1 + abs(z)));
krToKrcPrime = @(z) exp(0.25j * pi ./ (1 + abs(z))) ...
    .* (abs(z) + abs(z).^3 + (2 - 0.25j*pi)*z.^2) ...
    ./ abs(z) ./ (1 + abs(z)).^2;

krcToKr = @(z) z .* exp(-0.25j * pi ./ (1 + abs(z)));
krcToKrPrime = @(z) exp(-0.25j * pi ./ (1 + abs(z))) ...
    .* (abs(z) + abs(z).^3 + (2 + 0.25j*pi)*z.^2) ...
    ./ abs(z) ./ (1 + abs(z)).^2;


%% Integral
I1 = integral(f, 0, inf);
I2 = integral(@(x) f(krToKrc(x)) .* krToKrcPrime(x), 0, inf);
I3 = integral(@(x) f(krcToKr(x)) .* krcToKrPrime(x), 0, inf);

%% Plotting
rPlot(:, 1) = linspace(-2, 2, 1000);
iPlot(1, :) = linspace(-2, 2, 1000);

figure;
showImage(rPlot, iPlot, krToKrc(rPlot + 1j*iPlot), DisplayFormat="Real");
axis normal;

figure;
showImage(rPlot, iPlot, krToKrc(rPlot + 1j*iPlot), DisplayFormat="Imag");
axis normal;

figure;
showImage(rPlot, iPlot, krToKrcPrime(rPlot + 1j*iPlot), DisplayFormat="Real");
axis normal;

figure;
showImage(rPlot, iPlot, krToKrcPrime(rPlot + 1j*iPlot), DisplayFormat="Imag");
axis normal;




