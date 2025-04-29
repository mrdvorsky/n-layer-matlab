clc;
clear;
close all;

%% Inputs
NL = nLayerCircularTE(2, waveguideR=29/128*25.4);

f(:, 1) = linspace(32, 40, 1001);

er = 4;
thk = 2;

%% Get Cutoffs
kc1 = NL.modeCutoffs(1);
kc2 = NL.modeCutoffs(2);

k0 = 2*pi * f ./ 299.792468;

k1 = sqrt(k0.^2 .* er - kc1.^2);
k2 = sqrt(k0.^2 .* er - kc2.^2);

%% Plotting
figure;
plot(f, rad2deg(2*thk*k1), "", LineWidth=1.5);
hold on;
plot(f, rad2deg(2*thk*k2), "", LineWidth=1.5);
plot(f, rad2deg(2*thk*k0), "", LineWidth=1.5);
grid on;
ylim([0, inf]);







