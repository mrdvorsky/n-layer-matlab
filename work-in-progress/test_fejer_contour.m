clc;
clear;
close all;

%% Inputs
w = @(x) sinc(x).^2 ./ (1 + x.^2);

k = 5;
L = 2;

%% Test
f1 = @(x) w(x) .* exp(2j*k * acot(sqrt((x)./L))) ...
               ./ sin(2   * acot(sqrt((x)./L)));
M_k1 = integral(f1, 0, inf)

% f2 = @(x) (2*L) * w(L * cot(0.5*x).^2) ...
%                .* sin(k*x) ./ (1 - cos(x)).^2;
f2 = @(x) (2*L) * w(L * cot(0.5*x).^2) ...
               .* exp(1j*k*x) ./ (1 - cos(x)).^2;
M_k2 = integral(f2, 0, pi)


M_k2_2 = integral(f2, 0.0000001j, 3j) + integral(f2, pi + 3j, pi)


%% Contour Map
x1(:, 1) = linspace(-10, 20, 1000);
y1(1, :) = linspace(-2, 2, 1000);

v1 = f1(x1 + 1j*y1);

figure;
showImage(x1, y1, v1, DisplayFormat="Phase");
grid on;
axis normal;

figure;
showImage(x1, y1, v1, DisplayFormat="dB");
grid on;
axis normal;


x2(:, 1) = linspace(-0.5, 0.5, 1001);
y2(1, :) = linspace(-0.5, 0.5, 1001);

v2 = f2(x2 + 1j*y2);

figure;
showImage(x2, y2, v2, DisplayFormat="Phase");
grid on;
axis normal;

figure;
showImage(x2, y2, v2, DisplayFormat="dB");
grid on;
axis normal;
clim([-150, 150]);




