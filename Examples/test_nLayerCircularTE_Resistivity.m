% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
wgR = 5.8;

f = linspace(32, 40, 201).';

sheetRes = [10, 50, 100, 500, 1000];
sheetRes = linspace(0.8, 1.2, 5) .* 10;
sheetThk = 0.1;
sheetSigma = 1 ./ (sheetRes .* sheetThk);

numModes = 3;

%% Create nLayer Object
NL = nLayerCircularTE(numModes, waveguideR=wgR);

%% Calculate structure 4
[sheetEr] = NL.changeStructureConductivity(f, 1, 1, inf, sheetSigma(1));

er = [9 - 0.05j + 0*sheetEr, sheetEr];
thk = [0.7, sheetThk];

NL.printStructure(er(1, :), [], thk, ...
    Title=sprintf("Sheet Resistivity = %g Ohms/sq", sheetRes(1)));

tic;
gam = zeros(length(f), length(sheetRes));
for ii = 1:size(gam, 2)
    [er(:, 2)] = NL.changeStructureConductivity(f, 1, 1, inf, sheetSigma(ii));
    gam(:, ii) = NL.calculate(f, er, [], thk);
end
fprintf("Resonant case various conductivities: ");
toc;

% Plot
figure;
plot(gam, "-", Linewidth=1.5);
hold on;
zplane([]);
legend(compose("Sheet Resistivity = %g Ohms/sq", sheetRes), Location="southwest");
grid on;
