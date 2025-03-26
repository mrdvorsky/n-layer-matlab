clc;
clear;
close all;

%% Inputs
orderNum = 4;

xi1(:, 1) = linspace(0, 1, 21);
xi2 = linspace(0, 1, 21) .* (1 - xi1);

xt(1, 1, :) = [0, 1, 0.5];
yt(1, 1, :) = [0, 0, 1];

%% Get Shape Functions
shapeFun = fem_generateShapeFunctions(orderNum);

%% Get Shape Function Values
[val] = fem_evaluateShapeFunctions(shapeFun, xi1, xi2);

%% Get Global Coordinates
xi123 = cat(3, xi1 + 0*xi2, xi2, 1 - xi1 - xi2);
x = sum(xi123 .* xt, 3);
y = sum(xi123 .* yt, 3);

%% Get Coordinates of Zeros
[xs, ys] = fem_getGlobalElementCoordinates(orderNum, xt, yt, [1, 2, 3]);

%% Plot
for ii = 1:size(val, 3)
    figure;
    surf(x, y, val(:, :, ii), EdgeColor="none");
    shading interp;
    hold on;
    plot(xs, ys, ".w", MarkerSize=10);
    axis image;
    colormap colormapPlusMinus;
    clim(max(abs(clim)) * [-1, 1]);
    view(0, 90);
end



