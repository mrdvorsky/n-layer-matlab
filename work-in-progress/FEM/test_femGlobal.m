clc;
clear;
close all;

%% Inputs
orderNum = 4;

x(:, 1) = [-1, 1, 0, 0];
y(:, 1) = [0, 0, 1, -1];
tris = [3, 2, 1; 1, 4, 2];

%% Global
[x, y, tris] = fem_getGlobalElementCoordinates(x, y, tris, orderNum);

%% Plot
figure;
triplot(tris(:, 1:3), x, y);
hold on;
plot(x(tris(1, :)), y(tris(1, :)), "o");
plot(x(tris(2, :)), y(tris(2, :)), "x");
axis image;
xlim([-1.1, 1.1]);
ylim([-1.1, 1.1]);

for ii = 1:numel(x)
    text(x(ii), y(ii), num2str(ii), ...
        HorizontalAlignment="center", FontSize=22);
end


