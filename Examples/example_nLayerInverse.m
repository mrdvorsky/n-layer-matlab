% This example file shows how to use nLayerViewer to look at the outputs
% of various nLayerForward calculators.
%
% The basic usage of nLayerViewer is as follows (for single calculator):
%   NL = nLayerRectangular(...);
%   nLayerViewer();
%
% In the above example, maxM and maxN are the maximum mode m and n indices
% for any considered TEmn and TMmn modes. Typically, they are set to 3 and
% 2, respectively. Also, f is a column vector of frequencies, and er, ur,
% and thk are row vectors of complex permittivity and permeability and
% thickness for each layer. Optionally, er and ur can be matrices where
% each row corresponds to a particular frequency.
%
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
NL.printStructure(er, ur, thk);
gamActual1 = NL.calculate(f, er, ur, thk);
gamMeas1 = gamActual1 + (sqrt(0.5) .* noiseStd) ...
    .* complex(randn(size(f)), randn(size(f)));

%% Solve for Structure
NLsolver = nLayerInverse(2, verbosity=1);
NLsolver.setLayersToSolve(Erp=[2], Erpp=[2], Urp=[2], Urpp=[2], Thk=[]);
NLsolver.setRanges(ErpMin=[1, 0.01]);
NLsolver.setInitialValues(Er=er, Ur=ur, Thk=thk);
NLsolver.useGlobalOptimizer = false;

NLsolver.printStructureParameters(ShowLimits=true, Title="Case 1: Input");

tic;
[Params, Gamma, Uncert] = NLsolver.solveStructure(NL, f, gamMeas1);
toc;

NLsolver.printStructureParameters(Params, Uncert, Title="Case 1: Output");

NLsolver.computeParameterUncertainty(NL, f, NoiseStd=noiseStd);

%% Plot
figure;
nLayerViewer(Params.er, Params.thk, NL, f);
hold on;
plot(gamMeas1, "", Linewidth=1.5);
legend("Fit", "Measured");

