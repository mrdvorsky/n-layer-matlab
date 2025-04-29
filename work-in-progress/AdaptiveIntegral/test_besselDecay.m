clc;
clear;
close all;

%% Inputs
r = linspace(0.001, 30, 1001);

n = 0;
kc = besselj_zeros(n, 2);
% kc = kc(end);

c_start = 1;
rc = 0.4*max(r - c_start, 0);

%% Calculate
val1 = besselj(n, kc(end - 1) .* r);

val2 = 0.5 * (besselh(n, 1, kc(end) .* r + 1j*rc) + besselh(n, 2, kc(end) .* r - 1j*rc));
val3 = 0.5 * (besselh(n, 1, kc(end - 1) .* r + 1j*rc) + besselh(n, 2, kc(end - 1) .* r - 1j*rc));

%% Plotting
figure;
plot(r, real(val1), "", LineWidth=1.5);
hold on;
plot(r, imag(val1), "", LineWidth=1.5);
grid on;

figure;
plot(r, real(val2), "", LineWidth=1.5);
hold on;
plot(r, imag(val2), "", LineWidth=1.5);
grid on;

figure;
plot(r, real(val3), "", LineWidth=1.5);
hold on;
plot(r, imag(val3), "", LineWidth=1.5);
grid on;



