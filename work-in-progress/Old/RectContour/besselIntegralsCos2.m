clc;
clear;
close all;

%% Inputs
k(:, 1) = linspace(0.001, 20, 1001);

kp = 0.8;

%% Integral
for kk = 1:numel(k)
    kk
    jInt(kk) = integral(@(x) cos(x).^2 .* exp(-1j .* k(kk) .* cos(x - kp)), ...
        0, 2*pi);
end

%% Plotting
figure;
plot(k, real(jInt), "", LineWidth=1.5);
hold on;
plot(k, imag(jInt), "", LineWidth=1.5);
plot(k, pi*besselj(0, k) - cos(2*kp) *pi*besselj(2, k), ":", LineWidth=1.5);
grid on;








