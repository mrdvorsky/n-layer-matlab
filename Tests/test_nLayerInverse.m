% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Create Measurement Data
f = linspace(26.5, 40, 21).';
er = [1 - 0.0001j, 4 - 0.05j];
thk = [0.5, 0.5];
noiseStd = 0.01;

NL = nLayerRectangular(3, 2, waveguideBand="ka");
NL.printStructure(er, [], thk);
gamActual = NL.calculate(f, er, [], thk);
gamMeas1 = gamActual + (sqrt(0.5) .* noiseStd) ...
    .* complex(randn(size(f)), randn(size(f)));

%% Solve for Structure
NLsolver = nLayerInverse(2, verbosity=0);
NLsolver.setLayersToSolve(Erp=[2], Erpp=[], Thk=[1]);
NLsolver.setInitialValues(Er=er, Thk=thk);
NLsolver.useGlobalOptimizer = false;

NLsolver.rangeMax_thk(1) = 0.5;

NLsolver.printStructureParameters(ShowLimits=true, Title="Case 1: Input");

tic;
[er, ur, thk, gam] = NLsolver.solveStructure(NL, f, gamMeas1);
toc;

NLsolver.printStructureParameters(er, ur, thk, Title="Case 1: Output");

% NLsolver.computeParameterUncertainty(NL, f, NoiseStd=noiseStd);

%% Plot
% figure;
% nLayerViewer(er, thk, NL1, f1);
% hold on;
% plot(gamMeas1, "", Linewidth=1.5);
% legend("Fit", "Measured");

