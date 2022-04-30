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

%% Example 1: Ka-band rectangular waveguide
f1 = linspace(26.5, 40, 41).';
er1 = [1 - 0.01j, 4 - 0.05j];
thk1 = [0.5, 0.5];
noiseStd = 0.01;

NL1 = nLayerRectangular(3, 2, Band="ka");
gamActual1 = NL1.calculate(f1, er1, [], thk1);
gamMeas1 = gamActual1 + noiseStd.*randn(size(f1));

NLsolver = nLayerInverse(2, Verbosity=true);
NLsolver.useGlobalOptimizer = false;
NLsolver.setLayersToSolve(ErLayers=[2], ErpLayers=[2], ThkLayers=[1]);
NLsolver.setInitialGuesses(ErGuess=[1, 1], ErpGuess=[0.01, 0.05], ...
    ThkGuess=[0.01, 0.05]);

NLsolver.printStructureSetup(ShowLimits=false);

% tic;
% [er, ur, thk] = NLsolver.solveStructure(NL1, f1, gamMeas1);
% toc;





