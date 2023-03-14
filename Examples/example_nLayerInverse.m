% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
f = linspace(26.5, 40, 21).';
er = [1, 2 - 0.05j];
ur = [1, 0.5 - 0.01j];
thk = [0.5, 0.5];

noiseStd = 0.01;

%% Create Measurement Data
NL = nLayerFilledRectangular(1, 0, waveguideBand="ka", checkStructureValues=false);
gamActual = NL.calculate(f, er, ur, thk);
gamMeas = gamActual + (sqrt(0.5) .* noiseStd) ...
    .* complex(randn(size(f)), randn(size(f)));

%% Solve for Structure
NLsolver = nLayerInverse(2, verbosity=1);
NLsolver.setLayersToSolve(Erp=[2], Erpp=[2], Urp=[2], Urpp=[2], Thk=[]);
NLsolver.setRanges(UrpMin=[1, 0.01]);
NLsolver.setInitialValues(Er=er, Ur=ur, Thk=thk);
NLsolver.useGlobalOptimizer = false;

NLsolver.printStructureParameters(ShowLimits=true, Title="Input");

tic;
[Params, Gamma, Uncert] = NLsolver.solveStructure(NL, f, gamMeas);
toc;

% Uncert = NLsolver.computeParameterUncertainty(NL, f, NoiseStd=noiseStd);

NLsolver.printStructureParameters(Params, Uncert, Title="Output");

%% Plot
figure;
nLayerViewer(Params.er, Params.ur, Params.thk, NL, f);
hold on;
h = plot(gamMeas(:, :), ":", Linewidth=1.5);
set(h, {'DisplayName'}, cellstr(compose("%s (Meas)", NL.getOutputLabels().')));

