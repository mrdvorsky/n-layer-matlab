% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs


%% Create Measurement Data
f = linspace(26.5, 40, 21).';
er = [1, 0.5 - 0.05j];
ur = [1, 1 - 0.01j];
thk = [0.5, 0.5];

noiseStd = 0.01;

NL = nLayerRectangular(3, 2, waveguideBand="ka", checkStructureValues=false);
gamActual1 = NL.calculate(f, er, ur, thk);
gamMeas1 = gamActual1 + (sqrt(0.5) .* noiseStd) ...
    .* complex(randn(size(f)), randn(size(f)));

%% Solve for Structure
NLsolver = nLayerInverse(2, verbosity=1);
NLsolver.setLayersToSolve(Erp=[2], Erpp=[2], Urp=[2], Urpp=[2], Thk=[]);
NLsolver.setRanges(ErpMin=[1, 0.01]);
NLsolver.setInitialValues(Er=er, Ur=ur, Thk=thk);
NLsolver.useGlobalOptimizer = false;

NLsolver.printStructureParameters(ShowLimits=true, Title="Input");

tic;
[Params, Gamma, Uncert] = NLsolver.solveStructure(NL, f, gamMeas1(:, :));
toc;

Uncert = NLsolver.computeParameterUncertainty(NL, f, NoiseStd=noiseStd);
NLsolver.printStructureParameters(Params, Uncert, Title="Output");

%% Plot
figure;
nLayerViewer(Params.er, Params.ur, Params.thk, NL, f);
hold on;
h = plot(gamMeas1(:, :), ":", Linewidth=1.5);
set(h, {'DisplayName'}, cellstr(compose("%s (Meas)", NL.getOutputLabels().')));

