% This example file shows how to use nLayerCircularTM to calculate S11 of
% a circular waveguide looking into a multilayer structure.
%
% The basic usage of nLayerCircularTM is as follows (example for ka-band):
%   NL = nLayerRectangular(numModes, waveguideR=4.5);
%   S11 = NL.calculate(f, er, ur, thk);
%
% In the above example, numModes is the number of TM0n modes to consider.
% Typically, 3 modes are sufficient. Also, f is a column vector of
% frequencies, and er, ur, and thk are row vectors of complex permittivity
% and permeability and thickness for each layer. Optionally, er and ur can
% be matrices where each row corresponds to a particular frequency.
%
% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
wgR = 4.5;

% More frequency samples for the resonant case
f = linspace(26.5, 40, 51).';
fRes = linspace(26.5, 40, 2001).';

% Number of TM0n modes to consider
numModes = 3;

% Tolerance for convergence (-60 dB)
convergenceTol = 1e-3;

%% Create nLayer Object
NL = nLayerCircularTM(numModes, waveguideR=wgR);
% NL.interpolationPoints_kRho = 2.^16;
% NL.recomputeInterpolants();

% Optionally, specify convergence tolerance and verbosity
NL.convergenceAbsTol = convergenceTol;
NL.verbosity = 1;

%% Calculate structure 1
% Two layer conductor-backed structure (very low loss).
er1 = [1 - 0.0001j, 2 - 0.0001j];
ur1 = [];
thk1 = [1, 5];

NL.printStructure(er1, ur1, thk1, Title="Case 1");
tic;
gam1 = NL.calculate(f, er1, ur1, thk1);
fprintf("Two-layer low loss conductor backed: ");
toc;

% Plot
figure;
plot(gam1, ".-", Linewidth=1);
hold on;
title("Case 1");
zplane([]);
grid on;

%% Calculate structure 2
% Two layer infinte half-space structure
er2 = [1 - 0.0001j, 2 - 0.0001j];
ur2 = [];
thk2 = [1, inf];

NL.printStructure(er2, ur2, thk2, Title="Case 2");
tic;
gam2 = NL.calculate(f, er2, ur2, thk2);
fprintf("Two-layer low loss half-space: ");
toc;

% Plot
figure;
plot(gam2, ".-", Linewidth=1);
hold on;
title("Case 2");
zplane([]);
grid on;

%% Calculate structure 3
% One layer infinte half-space structure with permeability
er3 = [4 - 0.0001j];
ur3 = [1 - 0.0j];
thk3 = [inf];

NL.printStructure(er3, ur3, thk3, Title="Case 3");
tic;
gam3 = NL.calculate(f, er3, ur3, thk3);
fprintf("One-layer lossy half-space: ");
toc;

% Plot
figure;
plot(gam3, ".-", Linewidth=1);
hold on;
title("Case 3");
zplane([]);
grid on;

%% Calculate structure 4
% One-layer conductor backed resonant case
er4 = [4 - 1j];
ur4 = [];
thk4 = [10];

NL.printStructure(er4, ur4, thk4, Title="Case 4");
tic;
gam4 = NL.calculate(fRes, er4, ur4, thk4);
fprintf("One-layer conductor backed resonant: ");
toc;

% Plot
figure;
plot(gam4, ".-", Linewidth=1);
hold on;
title("Case 4");
zplane([]);
grid on;

%% Calculate structure 5
% Two-layer conductor backed frequency-variable structure
er5 = [linspace(1 - 0.0001j, 1 - 1j, length(f)).', ...
    linspace(2 - 0.0001j,4 - 4j, length(f)).'];
ur5 = [];
thk5 = [1, 5];

NL.printStructure(er5(1, :), ur5, thk5, Title="Case 5");
tic;
gam5 = NL.calculate(f, er5, ur5, thk5);
fprintf("One-layer conductor backed resonant: ");
toc;

% Plot
figure;
plot(gam5, ".-", Linewidth=1);
hold on;
title("Case 5");
zplane([]);
grid on;


