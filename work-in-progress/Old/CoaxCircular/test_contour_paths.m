clc;
clear;
close all;

%% Inputs
phi(:, 1) = linspace(0, 1, 102);
phi = phi(2:end-1, 1);

krP = phi;
kr = (1 - krP) ./ krP;

r = complex(kr, kr.^0.5 .* (kr < 1) + (kr >= 1) ./ 1);
r = complex(kr, kr.^2 ./ (kr.^2 + 1));
rP = 1 ./ (1 - r);

%% Calculate
% kr = L * (1 - krP) ./ krP;

%% Plotting
figure;
plot(complex(r), "o-", LineWidth=1.5);
hold on;
zplane([]);
grid on;
% xlim([-1, 10]);


figure;
plot(complex(rP), "o-", LineWidth=1.5);
hold on;
zplane([]);
grid on;





