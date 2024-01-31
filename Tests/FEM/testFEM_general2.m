clc;
clear;
close all;

%% Inputs
wgA = 1;
wgB = 0.5;

numX = 30;
numY = 25;

numModes = 10;

orderNum = 1;
orderNum = 2;
orderNum = 3;
orderNum = 4;
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

%% FEM Solution
[vals, cutoffs, x, y, tris, isBoundaryPEC, C, T] = fem_solveModes(x, y, tris, isBoundaryPEC, ...
    ModeType="TE", ModeCount=numModes, FemElementOrder=orderNum);
% db(cutoffs - pi*[1, 2, 2, hypot(1, 2), hypot(2, 2)].')

% [V, D] = eigs(C, T, numModes + 1, ...
%     "smallestabs", IsSymmetricDefinite=true);
% vals = V(:, 2:end);
% cutoffs = sqrt(diag(D(2:end, 2:end)));

%% Redo Triangulation for Plotting
trisPlot = delaunay(x, y);
trisEgdes = edges(triangulation(trisPlot, x, y));

%% Plotting
figure;
triplot(tris(:, 1:3), x, y);
hold on;
triplot(trisPlot, x, y, ":r");
plot(x(isBoundaryPEC), y(isBoundaryPEC), "ok", MarkerSize=7);
% plot(x(1 + numNonBoundary:end), y(1 + numNonBoundary:end), "o");
axis image;

for ii = 1:size(vals, 2)
    figure;
    trisurf(trisPlot(:, 1:3), x, y, 0*vals(:, ii), real(vals(:, ii)));
    colormap colormapPlusMinus;
    shading interp;
    hold on;
%     triplot(DT, x, y, 'Color', 'black');
    axis image;
    view(0, 90);
end

