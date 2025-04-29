clc;
clear;
close all;

%% Inputs
Nc = 50;
L = 1;
Lc = 1.0;

Nrho = 16000;
Nphi = 450;

%% Structure
er = 100 - 0.000j;
ur = 1;
thk = 2.5*20;
f = 30;
k0 = f * (2*pi/299.792458);

%% nLayer
wgA = 1;
wgB = 0.5;

%% Data
rPlot(:, 1) = 10 * linspace(-1, 1, 3001);
iPlot(1, :) = 10 * linspace(-1, 1, 3001);
cPlot = rPlot + 1j*iPlot;

gH = gamH(cPlot, k0, er, ur, thk);
gE = gamE(cPlot, k0, er, ur, thk);

%% Plotting
figure;
showImage(rPlot, iPlot, gH, DisplayFormat="dB");
clim([-50, 50]);

figure;
showImage(rPlot, iPlot, gE, DisplayFormat="dB");
clim([-50, 50]);






