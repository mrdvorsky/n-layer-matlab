clc;
clear;
close all;

%% Inputs
f = @(xs, ys) 1e1 * (1 + xs.^2 + 1.1*ys.^2)*0 + 1;

kx = 10.2;
ky = 10.1;

Tx(:, 1) = 0.01 * [0, 1, 0.7] + 10;
Ty(:, 1) = 0.01 * [0, 1, 0.5] + 100;

x(:, 1) = linspace(min(Tx) - 0.1, max(Tx) + 0.1, 2001);
y(1, :) = linspace(min(Ty) - 0.1, max(Ty) + 0.1, 2001);

%% Sample Function
f_cropped = @(xs, ys) f(xs, ys) .* inpolygon(xs + 0*ys, 0*xs + ys, Tx, Ty);
% f_sampled = f_cropped(x, y);

%% Integral Function
g = @(xs, ys) f(xs, ys) .* exp(1j .* (kx.*xs + ky.*ys));
g_cropped = @(xs, ys) f_cropped(xs, ys) .* exp(1j .* (kx.*xs + ky.*ys));

%% Midpoints of Triangle
Tx_mid = 0.5 * (Tx + circshift(Tx, 1));
Ty_mid = 0.5 * (Ty + circshift(Ty, 1));

%% Perform Integrals
% int1 = integral2(g_cropped, min(x), max(x), min(y), max(y))
int2 = TriIntegral(g_cropped, Tx, Ty)
% int3 = sum(g_cropped(x, y) .* abs(x(1) - x(2)) .* abs(y(1) - y(2)), "all")
int4 = sum(g(Tx_mid, Ty_mid)) .* polyarea(Tx, Ty) ./ 3

relErr = abs((int2 - int4) ./ int2)

%% Plotting
% figure;
% surf(x, y, f_sampled.', LineStyle="none");
% hold on;
% plot3(Tx, Ty, f(Tx, Ty), "kx", LineWidth=1.5);
% plot3(Tx_mid, Ty_mid, f(Tx_mid, Ty_mid), "ko", LineWidth=1.5);










