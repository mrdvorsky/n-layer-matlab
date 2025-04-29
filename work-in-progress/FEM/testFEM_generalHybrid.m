clc;
clear;
close all;

%% Inputs
wgA = 1;
wgB = 0.5;

numX = 31;
numY = 25;

numModes = 10;

orderNum = 1;
orderNum = 2;
% orderNum = 3;
% orderNum = 4;
% orderNum = 5;

%% Create Grid Points
[x, y] = ndgrid(wgA * linspace(-0.5, 0.5, numX), ...
    wgB * linspace(-0.5, 0.5, numY));
x = x(:);
y = y(:);

isBoundaryPEC = (x == min(x) | x == max(x) | y == min(y) | y == max(y));

x(~isBoundaryPEC) = x(~isBoundaryPEC) + 0.000 * (rand(size(x(~isBoundaryPEC))) - 0.5);
y(~isBoundaryPEC) = y(~isBoundaryPEC) + 0.000 * (rand(size(y(~isBoundaryPEC))) - 0.5);

% Perform Triangulation.
tris = delaunay(x, y);

er = zeros(size(tris, 1), 1);
er(mean(x(tris), 2) > 0) = 4;

%% FEM Solution
[valsE, valsH, cutoffs, x, y, tris, isBoundaryPEC, C, T] = fem_solveModesHybrid(x, y, tris, isBoundaryPEC, 0.5, ...
    ModeType="Hybrid", ModeCount=numModes, FemElementOrder=orderNum);
cutoffs

%% Redo Triangulation for Plotting
trisPlot = delaunay(x, y);
trisEgdes = edges(triangulation(trisPlot, x, y));

%% Plotting
figure;
triplot(tris(:, 1:3), x, y);
hold on;
triplot(trisPlot, x, y, ":r");
plot(x(isBoundaryPEC), y(isBoundaryPEC), "ok", MarkerSize=7);
axis image;

figure;
triplot(tris(:, 1:3), x, y);
hold on;
patch(x(tris(:, 1:3).'), y(tris(:, 1:3).'), real(er));
axis image;

for ii = 1:size(valsE, 2)
    figure;
    subplot(2, 1, 1);
    trisurf(trisPlot(:, 1:3), x, y, 0*valsE(:, ii), real(valsE(:, ii)));
    colormap colormapPlusMinus;
    shading interp;
    axis image;
    view(0, 90);

    subplot(2, 1, 2);
    trisurf(trisPlot(:, 1:3), x, y, 0*valsH(:, ii), real(valsH(:, ii)));
    colormap colormapPlusMinus;
    shading interp;
    axis image;
    view(0, 90);
end

