% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
f = linspace(26.5, 40, 201).';
er = [2 - 0.05j, 1];
ur = [0.5 - 0.01j, 1];
thk = [0.5, 0.5];

noiseStd = 0.01;

%% Create Measurement Data
NL = nLayerFilledRectangular(1, 0, waveguideBand="ka", checkStructureValues=false);
gamActual = NL.calculate(f, er, ur, thk);
gamMeas = gamActual + (sqrt(0.5) .* noiseStd) ...
    .* complex(randn(size(f)), randn(size(f)));
% gamMeas = gamActual + (sqrt(1) .* noiseStd) ...
%     .* exp(1j*2*pi*rand(1, 1));

% load tmp;
% gamMeas = gamActual + reshape(J(:, :, 4), [], 2, 2);

%% Solve for Structure
NLsolver = nLayerInverse(numel(thk), verbosity=1);
NLsolver.setLayersToSolve(Er=[1], Ur=[1], Thk=[1, 2]);
NLsolver.setRanges(UrpMin=[0.1]);
NLsolver.setInitialValues(Er=er, Ur=ur, Thk=thk);
% NLsolver.addThicknessConstraint("all", IsFixed=true);
NLsolver.addThicknessConstraint(1, IsFixed=true)
NLsolver.useGlobalOptimizer = false;
% NLsolver.localOptimizerOptions = optimoptions("fmincon");

NLsolver.printStructureParameters(ShowLimits=true, Title="Input");

tic;
[Params, Gamma, Uncert] = NLsolver.solveStructure(NL, f, gamMeas, NoiseStdMin=noiseStd);
toc;

Uncert = NLsolver.computeParameterUncertainty(NL, f, NoiseStd=noiseStd);

NLsolver.printStructureParameters(Params, Uncert, Title="Output");

%% Plot
figure;
nLayerViewer(Params.er, Params.ur, Params.thk, NL, f);
hold on;
h = plot(gamMeas(:, :), ":", Linewidth=1.5);
set(h, {'DisplayName'}, cellstr(compose("%s (Meas)", NL.getOutputLabels().')));

% J = reshape(full(Uncert.J(1:4*numel(f), :) + 1j*Uncert.J(4*numel(f)+1:end, :)), numel(f), 4, []);
% save tmp J;
% figure;
% plot(J(:, :, 1));


