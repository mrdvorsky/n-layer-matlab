% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
f = linspace(26.5, 40, 21).';
er = [2 - 0.05j, 1];
ur = [0.5 - 0.01j, 1];
thk = [0.5, 0.1];

noiseStd = 0.01;

%% Create Measurement Data
NL = nLayerFilledRectangular(1, 0, waveguideBand="ka", checkStructureValues=false);
gamActual = NL.calculate(f, er, ur, thk);
gamMeas = gamActual + (sqrt(0.5) .* noiseStd) ...
    .* complex(randn(size(f)), randn(size(f)));

%% Solve for Structure
NLsolver = nLayerInverse(numel(thk), verbosity=1);
NLsolver.setLayersToSolve(Erp=[1], Erpp=[1], Urp=[1], Urpp=[1], Thk=[]);
NLsolver.setRanges(UrpMin=[0.01, 1]);
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

