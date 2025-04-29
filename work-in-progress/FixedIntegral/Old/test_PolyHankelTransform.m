clc;
clear;
close all;

%% Inputs
p1 = [1, 0, 0];

kr(:, 1) = 50 * linspace(-1, 1, 1000);

ng = 1000;

a = [1/2];
b = [1];

%% Fixed Grid Method
[gr1(1, :), gw1(1, :)] = fejer2(ng, a, b);
spec1 = sum(gw1 .* besselj(0, kr .* gr1) .* gr1 .* polyval(p1, gr1), 2);

%% Analytical Method
[x0, c] = polyHankelTransform(a, b, ...
    [polyval(p1, linspace(a, b, numel(p1)))]);

spec2 = -(besselj(1, 1/3 * kr) - 3 * besselj(1, kr)) ./ (3 * kr);

%% Plot
plot(kr, real(spec1), "", LineWidth=1.5);
hold on;
grid on;
yL = ylim();

figure;
plot(kr, real(spec2), "", LineWidth=1.5);
hold on;
grid on;
ylim(yL);
 


