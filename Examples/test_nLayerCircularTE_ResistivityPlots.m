% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
% close all;

%% Inputs
wgR = 5.8;

f = linspace(32, 40, 201).';

sheetRes = logspace(log10(0.1), log10(1000), 101);
sheetThk = 0.1;
sheetSigma = 1 ./ (sheetRes .* sheetThk);

numModes = 3;

%% Create nLayer Object
NL = nLayerCircularTE(numModes, waveguideR=wgR);

%% Calculate structure 4
[sheetEr] = NL.changeStructureConductivity(f, 1, 1, inf, sheetSigma(1));

er = [9 - 0.05j + 0*sheetEr, sheetEr];
thk = [0.7, inf];

% er = [4 - 0.001j + 0*sheetEr, sheetEr];
% thk = [1.1, inf];

% er = [1 - 0.0j + 0*sheetEr, sheetEr];
% thk = [3.4, inf];

NL.printStructure(er(1, :), [], thk, ...
    Title=sprintf("Sheet Resistivity = %g Ohm/sq", sheetRes(1)));

tic;
gam = zeros(length(f), length(sheetRes));
for ii = 1:size(gam, 2)
    [er(:, 2)] = NL.changeStructureConductivity(f, 1, 1, inf, sheetSigma(ii));
    gam(:, ii) = NL.calculate(f, er, [], thk);
end
fprintf("Resonant case various conductivities: ");
toc;

%% Fitting
gamDist = min(abs(1 - gam), [], 1);

%% Plot
figure;
semilogx(sheetRes, gamDist, "-", LineWidth=1.5, ...
    DisplayName=sprintf("Applicator \\epsilon_r = %g", real(er(1))));
hold on;
xlabel("Sheet Resistivity, Ohm/sq");
ylabel("min(|1 - S_{11}|)");
grid on;

