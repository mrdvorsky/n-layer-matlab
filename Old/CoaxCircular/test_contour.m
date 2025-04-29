clc;
clear;
close all;

%% Inputs
krReal(:, 1) = 2*pi * linspace(-1, 1, 1001);
krImag(1, :) = 2*pi * linspace(-1, 1, 1001);

er = [1 - 0.01j];
thk = [3];

circR = 1;
intLength = 10;

intPoints = 1000000;

%% Calculate
[gam1, gam2] = nLayerOpenEnded.computeGamma0(krReal + 1j*krImag, 1, er, ones(size(er)), thk);

int1 = @(kr) nLayerOpenEnded.computeGamma0(kr, 1, er, ones(size(er)), thk) .* sinc(kr) .* exp(-imag(kr));

%% Integrals
x1 = (1 + 1j) * linspace(0, circR, intPoints);
x2 = (1j) * circR + linspace(circR, intLength, intPoints);
x3 = intLength + 1j * linspace(0, circR, intPoints);

y1 = mean(int1(x1)) * (1 + 1j) * circR;
y2 = mean(int1(x2)) * (intLength - circR);
y3 = mean(int1(x3)) * (-1j) * circR;

y = mean(int1(linspace(0, intLength, intPoints))) * intLength
y1 + y2 + y3

%% Plotting
figure;
showImage(krReal, krImag, gam1, DisplayFormat="dB");
colormap gray;
grid on;
% clim(max(abs(clim())) * [-1, 1]);
clim(max(abs(clim())) + [-80, 0]);

figure;
showImage(krReal, krImag, gam2, DisplayFormat="dB");
colormap gray;
grid on;
% clim(max(abs(clim())) * [-1, 1]);
clim(max(abs(clim())) + [-80, 0]);




