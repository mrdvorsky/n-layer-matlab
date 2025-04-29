% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
wgEr = 1 - 0.00j;
wgR = 5.8 ./ sqrt(real(wgEr));

fRes = linspace(32, 40, 201).';

numModes = 3;

%% Create nLayer Object
NL = nLayerCircularTE(numModes, waveguideR=wgR, waveguideEr=wgEr);

%% Calculate structure 4
er = [4 - 0.001j];
thk = [1.2];
sigma = 57e3 * flip([0.00001, 0.0001, 0.001, 0.01, 0.1, 1]);

NL.printStructure(er(1, :), [], thk, BackingConductivity=sigma(1), ...
    Title=sprintf("Waveguide er = %g", real(wgEr)));

tic;
gam = zeros(length(fRes), length(sigma));
for ii = 1:length(sigma)
    gam(:, ii) = NL.calculate(fRes, er, [], thk, ...
        BackingConductivity=sigma(ii));
end
fprintf("Resonant case various conductivities: ");
toc;

% Plot
figure;
plot(gam, "-", Linewidth=1.5);
hold on;
title("Case 4");
zplane([]);
legend(compose("\\sigma = %.1e S/m", 1000*abs(sigma)));
legend(compose("%g %%IACS", 100*abs(sigma ./ 57e3)));
% legend("Smooth", "RMS: 0.5\delta", "RMS: \delta", "RMS: >> \delta");
grid on;
