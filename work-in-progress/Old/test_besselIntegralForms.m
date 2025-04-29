clc;
clear;
close all;

%% Inputs
x(:, 1) = 5 * linspace(-1, 1, 1001);
y(1, :) = 5 * linspace(-1, 1, 1001);

k = 10;

%% Calc
f = @(phi) exp(1j .* k .* sin(phi));

% v = sin(x + 1j*y);
v = f(x + 1j*y);

%% Integral
I1 = 2*pi .* besselj(0, k)

I2 = integral(f, 0, 2 * pi)

I3 = integral(f, 0, 2j) + integral(f, 2*pi + 2j, 2*pi) ...
    + integral(f, 2j, pi) + integral(f, pi, 2*pi + 2j)


%% Plotting
figure;
showImage(x, y, v, DisplayFormat="Phase");
colormap hsv;

figure;
showImage(x, y, v, DisplayFormat="dB");
colormap jet;






