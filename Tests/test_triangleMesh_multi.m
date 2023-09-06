clc;
clear;
close all;

%% Inputs
f = @(xs, ys) 1e1 * (1 + xs + 0.1*ys)*0.001 + 1;

filename = "three";
filename = "gauss4x4";
filename = "gauss8x8";
filename = "toms706_37";

kRho(:, 1) = 10.^linspace(-0, 2, 11);
kRho = 100;
kPhi(1, :) = linspace(0, 90, 101);

Tx(:, 1) = 1 * [0, 1, 0] + 0;
Ty(:, 1) = 1 * [0, 0, 1] + 0;

%% Guassian Grid Points
% gx = [0.5, 0.5, 0].';
% gy = [0, 0.5, 0.5].';
% gw = [1, 1, 1].' ./ 3;

% Read from file.
xData = readmatrix(strcat(filename, "_x"));
wData = readmatrix(strcat(filename, "_w"));

gx = xData(:, 1);
gy = xData(:, 2);
gw = wData(:, end);

%% Perform Integrals
kx = kRho .* cosd(kPhi);
ky = kRho .* sind(kPhi);

relErr1 = 0*kx;
intVal = 0*kx;
for ii = 1:numel(kx)
    g = @(xs, ys) f(xs, ys) .* exp(1j .* (kx(ii).*xs + ky(ii).*ys));

    int1 = TriIntegral(g, Tx, Ty);
    int2 = sum(gw .* g(gx, gy)) .* polyarea(Tx, Ty);

    intVal(ii) = int1;
    relErr1(ii) = abs((int1 - int2) ./ int1);
end

%% Plotting
% figure;
% loglog(kRho, relErr1, "", LineWidth=1.5);
% hold on;
% xlabel("k_{\rho}");
% ylabel("Relative Error");
% title(filename);
% grid on;
% xlim([min(kRho), max(kRho)]);
% ylim([1e-10, 1]);

figure;
plot(kPhi, real(intVal), "", LineWidth=1.5);
hold on;
plot(kPhi, imag(intVal), "--", LineWidth=1.5);
xlabel("k_{\phi}");
ylabel("Val");
title(filename);
grid on;
% xlim([min(kRho), max(kRho)]);
% ylim([1e-10, 1]);













