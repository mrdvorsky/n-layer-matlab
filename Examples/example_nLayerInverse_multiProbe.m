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
% Example 1 
f1 = linspace(32, 40, 21).';
er1 = [1 - 0.0j, 4 - 0.05j];
thk1 = [0.5, 0.5];
noiseStd = 0.01;

%% Create Measurement Data
% Example 1: Rectangular and circularTE
NL1_1 = nLayerRectangular(3, 2, Band="ka");
gamMeas1_1 = NL1_1.calculate(f1, er1, [], thk1) + (sqrt(0.5) .* noiseStd) ...
    .* complex(randn(size(f1)), randn(size(f1)));

NL1_2 = nLayerCircularTE(3, R=5.8);
gamMeas1_2 = NL1_2.calculate(f1, er1, [], thk1) + (sqrt(0.5) .* noiseStd) ...
    .* complex(randn(size(f1)), randn(size(f1)));

%% Solve for Structure
% Example 1
NLsolver = nLayerInverse(2, Verbosity=1);
NLsolver.setInitialValues(ErValue=real(er1), ErpValue=-imag(er1), ThkValue=thk1);
NLsolver.thkInitialValue(1) = [0.5];    % Make this a "guess"
NLsolver.setLayersToSolve(ErLayers=[2], ErpLayers=[2], ThkLayers=[1]);

NLsolver.printStructureParameters(ShowLimits=true, Title="Case 1: Input");
[er, ur, thk, gam1_1, gam1_2] = NLsolver.solveStructure(...
    NL1_1, f1, gamMeas1_1, ...
    NL1_2, f1, gamMeas1_2);
NLsolver.printStructureParameters(er, ur, thk, Title="Case 1: Output");

% NLsolver.computeParameterUncertainty(NL1_1, f1, NoiseStd=noiseStd);
% NLsolver.computeParameterUncertainty(NL1_2, f1, NoiseStd=noiseStd);
% NLsolver.computeParameterUncertainty(NL1_1, f1, NL1_2, f1, NoiseStd=noiseStd);

%% Plot
figure;
nLayerViewer(er, thk, NL1_1, f1, NL1_2, f1);
hold on;
plot(gamMeas1_1, "", Linewidth=1.5);
plot(gamMeas1_2, "", Linewidth=1.5);
legend("Fit Rectangular", "Fit Circular", "Measured Rect", "Measured Circ");

