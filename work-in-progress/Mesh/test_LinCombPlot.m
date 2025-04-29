clc;
clear;
close all;

%% Inputs
f = @(xs, ys) 1 + 1*xs + 0*ys;

kx(:, 1) = 10 * linspace(-1, 1, 401);
ky(1, :) = 10 * linspace(-1, 1, 401);

krLine(:, 1) = linspace(1, 20, 51);
% kphiLine = 20;
kphiLine(:, 1) = linspace(0.1, 179.9, 400);

fitOrder = 3;

% Fixed Quadrature
nx = 1*128;
ny = 1*128;

% Coordinate Transformation
coordTransX = [2, 0.3, 0];        % x1 = c(1)*x + c(2)*y + c(3);
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

%% Perform Analytical Line Fit
for pp = 1:numel(kphiLine)
    pp
    specLine = nufftn(f(gx, gy) .* gw, [gx, gy]./(2*pi), ...
        krLine .* [cosd(kphiLine(pp)), sind(kphiLine(pp))]);

    kxLine = krLine .* cosd(kphiLine(pp));
    kyLine = krLine .* sind(kphiLine(pp));

    fitBase = exp(-1j * (kxLine .* Tx.' + kyLine .* Ty.'));
    fitAll = [];
    for ii = 1:fitOrder
        fitAll = [fitAll, fitBase ./ krLine.^ii];
    end

    fitCoeff(:, pp) = fitAll \ specLine;
    errRms(pp) = rms(fitAll * fitCoeff(:, pp) - specLine);
end

%% Plot Fit Lines
% figure;
% plot(krLine, real(specLine), "", LineWidth=1.5);
% hold on;
% plot(krLine, imag(specLine), "", LineWidth=1.5);
% plot(krLine, real(fitAll * fitCoeff), "x", LineWidth=1.5);
% plot(krLine, imag(fitAll * fitCoeff), "x", LineWidth=1.5);
% grid on;

%% Plot Angle Lines
figure;
plot(kphiLine, real(fitCoeff(4:6, :).'), "", LineWidth=1.5);
grid on;

% figure;
% plot(kphiLine, imag(fitCoeff.'), "", LineWidth=1.5);
% grid on;






















