clc;
clear;
close all;

%% Inputs
p1 = [1, -2, 1, 1, -1];
p2 = 1*[-1.5, 0, 0, 0, 0.5];

kx(:, 1) = 50 * linspace(-1, 1, 1000);

ng = 1000;

a = [-1; 1];
b = [2; 1.9];

%% Fixed Grid Method
[gx1, gw1] = fejer2(ng, a(1), b(1));
[gx2, gw2] = fejer2(ng, a(2), b(2));
spec1 = nufft(polyval(p1, gx1) .* gw1, gx1 ./ (2*pi), kx) ...
    + nufft(polyval(p2, gx2) .* gw2, gx2 ./ (2*pi), kx);

%% Analytical Method
[x0, c] = polyFourierTransform(a, b, ...
    [polyval(p1, linspace(a(1), b(1), numel(p1))); ...
    polyval(p2, linspace(a(2), b(2), numel(p2)))]);

spec2 = sum(nufft(c, x0 ./ (2*pi), kx, 1) ...
    .* (1./kx).^(1:size(c, 2)), 2);

%% Plot
plot(kx, real(spec1), "", LineWidth=1.5);
hold on;
plot(kx, imag(spec1), "", LineWidth=1.5);
grid on;
yL = ylim();

figure;
plot(kx, real(spec2), "", LineWidth=1.5);
hold on;
plot(kx, imag(spec2), "", LineWidth=1.5);
grid on;
ylim(yL);
 


