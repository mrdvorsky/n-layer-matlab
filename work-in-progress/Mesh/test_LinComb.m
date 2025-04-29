clc;
clear;
close all;

%% Inputs
f = @(xs, ys) 1 + 0*xs + 0*ys;

kx(:, 1) = 10 * linspace(-1, 1, 401);
ky(1, :) = 10 * linspace(-1, 1, 401);

krLine(:, 1) = linspace(1, 50, 501);
kphiLine = 1;

fitOrder = 3;

% Fixed Quadrature
nx = 2*128;
ny = 2*128;

% Coordinate Transformation
coordTransX = [1, 1, 0];        % x1 = c(1)*x + c(2)*y + c(3);
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

%% Line Sample
specLine = nufftn(f(gx, gy) .* gw, [gx, gy]./(2*pi), ...
    krLine .* [cosd(kphiLine), sind(kphiLine)]);

% specLine = (exp(-1j * kxLine).*kyLine./(kyLine - kxLine) ...
%     - exp(-1j * kyLine).*kxLine./(kyLine - kxLine) - 1) ...
%     ./ (kxLine.*kyLine);

%% Perform Analytical
kr = hypot(kx, ky);
spec2 = exp(-1j * kx) - exp(-1j * ky) ./ kr;

%% Perform Analytical Fit
kxLine = krLine .* cosd(kphiLine);
kyLine = krLine .* sind(kphiLine);

% fit1 = exp(-1j * Tx.' * kxLine);
% fit2 = exp(-1j * kyLine);
% fit3 = exp(-0j * kyLine);
fitBase = exp(-1j * (kxLine .* Tx.' + kyLine .* Ty.'));
fitAll = [];
for ii = 1:fitOrder
    fitAll = [fitAll, fitBase ./ krLine.^ii];
end

fitCoeff = fitAll \ specLine

errRms = rms(fitAll * fitCoeff - specLine)

%% Plot Grids
% figure;
% plot([Tx; Tx(1)], [Ty; Ty(1)], "", LineWidth=1.5);
% hold on;
% plot(gx, gy, "o", LineWidth=1.5);
% grid on;

%% Plot Lines
figure;
plot(krLine, real(specLine), "", LineWidth=1.5);
hold on;
plot(krLine, imag(specLine), "", LineWidth=1.5);
plot(krLine, real(fitAll * fitCoeff), "x", LineWidth=1.5);
plot(krLine, imag(fitAll * fitCoeff), "x", LineWidth=1.5);
grid on;

% figure;
% plot(krLine, real(fitAll), "", LineWidth=1.5);
% hold on;
% plot(krLine, imag(fitAll), "", LineWidth=1.5);
% grid on;

%% Plotting
% figure;
% showImage(kx, ky, spec1 .* hypot(kx, ky), DisplayFormat="Magnitude");
% xlabel("k_x");
% ylabel("k_y");

% figure;
% showImage(kx, ky, spec2, DisplayFormat="Magnitude");
% xlabel("k_x");
% ylabel("k_y");















