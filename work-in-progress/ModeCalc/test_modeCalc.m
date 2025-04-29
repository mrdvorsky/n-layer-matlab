clc;
clear;
close all;

%% Inputs
wgA = 7.112;
wgB = 3.556;

Nx = 10;
Ny = 10;

kr(1, 1, :) = linspace(0, 2, 200);

kc11 = hypot(pi./wgA, pi./wgB);

%% Quad Points
[xg(1, :), xw(1, :)] = fejer2(Nx, 0, wgA);
[yg(1, :), yw(1, :)] = fejer2(Ny, 0, wgB);

xAll(1, :) = [xg, 0*yg + wgA, flip(xg), 0*yg] - 0.5*wgA;
yAll(1, :) = [0*xg, yg, 0*xg + wgB, flip(yg)] - 0.5*wgB;

%% Spectrum Points
theta(:, 1) = linspace(0, 2*pi, 8*(Nx + Ny) + 1);
theta = theta(1:end - 1);

A = exp(-1j * kr .* (cos(theta).*xAll + sin(theta).*yAll));
S = pagesvd(A);

%% Plotting
figure;
plots(kr, S(:, 1, :), "", LineWidth=1.5);
grid on;




