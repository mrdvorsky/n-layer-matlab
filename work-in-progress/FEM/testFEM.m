clc;
clear;
close all;

%% Inputs
wgA = 1;
wgB = 0.5;

numX = 5;
numY = 5;

numModes = 12;

%% Create Grid Points
[x, y] = ndgrid(wgA * linspace(-0.5, 0.5, numX), ...
    wgB * linspace(-0.5, 0.5, numY));
x = x(:);
y = y(:);

% Sort x and y so edges are last.
isBoundaryPEC = (x == min(x) | x == max(x) | y == min(y) | y == max(y));

% Perform Triangulation.
DT = delaunay(x, y);

%% FEM Solution
[vals, cutoffs] = solveModesFEM(x, y, DT, isBoundaryPEC, ModeType="TE");

%% Plotting
figure;
triplot(DT, x, y);
hold on;
plot(x(1 + numNonBoundary:end), y(1 + numNonBoundary:end), "o");
axis image;

for ii = 1:size(vals, 2)
    figure;
    trisurf(DT, x, y, 0*vals(:, ii), vals(:, ii));
    colormap jet;
    shading interp;
    hold on;
%     triplot(DT, x, y, 'Color', 'black');
    axis image;
    view(0, 90);
end

