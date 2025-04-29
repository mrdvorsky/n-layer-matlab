clc;
clear;
close all;

%% Inputs
krMax = 200;

Nkr = 4000;
Nxy = 1000;

%% Waveguide
NL = nLayerRectangular(1, 0, waveguideBand="X");
specFun = NL.modeStructs(1).EySpec;

%% Sample
wgA = NL.waveguideA;
wgB = NL.waveguideB;

x(:, 1) = wgA * linspace(-1, 1, Nxy);
y(1, :) = wgB * linspace(-1, 1, Nxy);

%% Calculate Spectrum
kx(:, 1) = linspace(-krMax, krMax, Nkr);
ky(1, :) = linspace(-krMax, krMax, Nkr);

specEy = specFun(kx, ky);

xs = x + 0*y;
ys = 0*x + y;
scale = abs(x(2) - x(1)) * abs(y(2) - y(1));
Ey = reshape(scale ...
    .* nufftn(specEy, {kx, ky}, [xs(:), ys(:)] ./ (2*pi)), ...
    [numel(x), numel(y)]);

%% Theory
EyTheory = cos(x * pi ./ wgA) ...
    .* (abs(x) <= 0.5*wgA) .* (abs(y) <= 0.5*wgB);

%% Plotting
figure;
showImage(x, y, Ey, DisplayFormat="Magnitude");

figure;
showImage(x, y, EyTheory, DisplayFormat="Magnitude");





