clc;
clear;
close all;

%% Inputs
f = @(xs, ys) 1 + 0*xs + 0*ys;

kx(:, 1) = 30 * linspace(-1, 1, 401);
ky(1, :) = 30 * linspace(-1, 1, 401);

% Fixed Quadrature
nx = 128;
ny = 128;

% Coordinate Transformation
coordTransX = [1, 0, 0];        % x1 = c(1)*x + c(2)*y + c(3);
coordTransY = [0, 1, 0];        % y1 = c(1)*x + c(2)*y + c(3);

%% Generate Triangle Coordinates
Tx(:, 1) = coordTransX(3) + [0, coordTransX(1), coordTransX(2)];
Ty(:, 1) = coordTransY(3) + [0, coordTransY(1), coordTransY(2)];

%% Guassian Grid Points
[nodesX(:, 1), weightsX(:, 1)] = fejer2(nx, 0, 1);
[nodesY(1, :), weightsY(1, :)] = fejer2(ny, 0, 1);

gxUnit = nodesX + 0*nodesY;
gyUnit = nodesY .* (1 - nodesX);
gwUnit = weightsY .* weightsX .* (1 - nodesX);

gx = coordTransX(1).*gxUnit(:) + coordTransX(2).*gyUnit(:) + coordTransX(3);
gy = coordTransY(1).*gxUnit(:) + coordTransY(2).*gyUnit(:) + coordTransY(3);
gw = gwUnit(:) .* polyarea(Tx, Ty) .* 2;

%% Perform Fixed Quadrature
spec1 = reshape(...
    nufftn(f(gx, gy) .* gw, [gx, gy]./(2*pi), {kx, ky}), ...
    numel(kx), []);

%% Perform Analytical
% spec2 = exp(-1j * kx) - exp(-1j * ky) ./ kr;

kx = kx + 0.001;
ky = ky + 0.002;
kr = hypot(kx, ky);
spec2 = (exp(-1j * kx).*ky./(ky - kx) - exp(-1j * ky).*kx./(ky - kx) - 1) ...
    ./ (kx.*ky);

spec2 = exp(-1j * kx);

%% Plot Grids
% figure;
% plot([Tx; Tx(1)], [Ty; Ty(1)], "", LineWidth=1.5);
% hold on;
% plot(gx, gy, "o", LineWidth=1.5);
% grid on;

%% Plotting
figure;
showImage(kx, ky, spec1, DisplayFormat="Magnitude");
xlabel("k_x");
ylabel("k_y");

figure;
showImage(kx, ky, spec2, DisplayFormat="Magnitude");
xlabel("k_x");
ylabel("k_y");















