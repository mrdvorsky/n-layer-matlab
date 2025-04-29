% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
f = linspace(26.5, 40, 21).';
er = [1, 2 - 0.05j];
ur = [1, 0.5 - 0.01j];
thk1 = [0.5, 0.5];      % MUT 0.5 mm
thk2 = [0.5, 1.5];      % MUT 1.5 mm

noiseStd = 0.01;

%% Create Measurement Data
NL = nLayerFilledRectangular(1, 0, waveguideBand="ka", checkStructureValues=false);
gamActual1 = NL.calculate(f, er, ur, thk1);
gamActual2 = NL.calculate(f, er, ur, thk2);
gamMeas1 = gamActual1 + (sqrt(0.5) .* noiseStd) ...
    .* complex(randn(size(f)), randn(size(f)));
gamMeas2 = gamActual2 + (sqrt(0.5) .* noiseStd) ...
    .* complex(randn(size(f)), randn(size(f)));

%% Solve for Structure
NLsolver1 = nLayerInverse(2, verbosity=1);
NLsolver1.setLayersToSolve(Erp=[2], Erpp=[2], Urp=[2], Urpp=[2], Thk=[]);
NLsolver1.setRanges(UrpMin=[1, 0.01]);
NLsolver1.setInitialValues(Er=er, Ur=ur, Thk=thk1);
NLsolver1.useGlobalOptimizer = false;

NLsolver2 = copy(NLsolver1);
NLsolver2.initialValue_thk(:) = thk2;

NLsolver1.printStructureParameters(ShowLimits=true, Title="Input 1");
NLsolver2.printStructureParameters(ShowLimits=true, Title="Input 2");

tic;
[Params, Gamma, Uncert] = nLayerInverse.solveStructureMultiple(...
    NLsolver1, NL, f, gamMeas1, ...
    NLsolver2, NL, f, gamMeas2);
toc;

% Uncert = nLayerInverse.computeParameterUncertaintyMultiple(...
%     NLsolver1, NL, f, ...
%     NLsolver2, NL, f, ...
%     NoiseStd=noiseStd);

NLsolver1.printStructureParameters(Params{1}, Uncert{1}, Title="Output 1");
NLsolver2.printStructureParameters(Params{2}, Uncert{2}, Title="Output 2");

%% Plot
figure;
nLayerViewer(Params{1}.er, Params{1}.ur, Params{1}.thk, NL, f);
hold on;
h = plot(gamMeas1(:, :), ":", Linewidth=1.5);
set(h, {'DisplayName'}, cellstr(compose("%s (Meas)", NL.getOutputLabels().')));

figure;
nLayerViewer(Params{2}.er, Params{2}.ur, Params{2}.thk, NL, f);
hold on;
h = plot(gamMeas2(:, :), ":", Linewidth=1.5);
set(h, {'DisplayName'}, cellstr(compose("%s (Meas)", NL.getOutputLabels().')));

