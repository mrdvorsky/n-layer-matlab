clc;
clear;
close all;

%% Inputs
frac = 1;
x(:, 1) = linspace(0, 1, 100001);

order = 4;

%% Test Fit
y = exp(1j .* (2*pi) * frac * x) .* (1 + x + x.^2);
p = polyfit(x, y, order);

yFit = polyval(p, x);

%% Numerical
int1 = trapz(x, y)
int2 = trapz(x, yFit)

relErr = abs((int1 - int2) ./ int1)

%% Plot
figure;
plot(x, rad2deg(angle(y)), "", LineWidth=1.5);
hold on;
plot(x, rad2deg(angle(yFit)), "", LineWidth=1.5);




