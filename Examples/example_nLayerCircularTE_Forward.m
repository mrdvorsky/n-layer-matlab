% This example file shows how to use nLayerCircularTE to calculate S11 of
% a circular waveguide looking into a multilayer structure.
%
% This example mostly shows 1-layer cases and finite conductivity cases.
% For multilayer cases and frequency dependent cases, see the example for
% nLayerRectangular. These examples can be applied directly to
% nLayerCircularTE simply by creating the NL object with nLayerCircularTE
% instead.
%
% The basic usage of nLayerCircularTE is as follows (example for ka-band):
%   NL = nLayerRectangular(numModes, waveguideR=5.8);
%   S11 = NL.calculate(f, er, ur, thk);
%
% In the above example, numModes is the number of TE0n modes to consider.
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
% More frequency samples for the resonant case
f = linspace(32, 40, 51).';
fRes = linspace(32, 40, 2001).';

% Number of TE0n modes to consider
numModes = 3;

% Tolerance for convergence (-60 dB)
convergenceTol = 1e-3;

%% Create nLayer Object
% Creating the object initially takes a short amount of time (~0.01
% seconds) to precompute integral stuff. You can either create this once at
% the beginning of your program, or save the nLayer object to a file and
% load it (second line here). See documentation for nLayerCircularTE for
% more details.

NL = nLayerCircular(0, numModes, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
% NL = load("nLayer_circularKaBand").NL;

% Optionally, specify convergence tolerance and verbosity
% NL.convergenceAbsTol = convergenceTol;
% NL.verbosity = 1;

%% Calculate structure 1
% One layer conductor-backed structure (different losses).
% Layer 1: er = 4 - 0.0001j, ur = 1, thk = 2 mm
er1 = 4 - 1j.*10.^linspace(-0, -3, 7).';
ur1 = [];
thk1 = [1.2];

nLayer.printStructure(er1(1, :), ur1, thk1, Title="Case 1");
tic;
gam1 = zeros(length(f), length(er1));
for ii = 1:length(er1)
    gam1(:, ii) = NL.calculate(f, er1(ii), ur1, thk1);
end
fprintf("Various loss factors: ");
toc;

% Plot
figure;
plot(gam1, ".-", Linewidth=1);
hold on;
title("Case 1");
zplane([]);
legend(compose("\\epsilon'' = %.1e", abs(imag(er1))));
grid on;

%% Calculate structure 2
% One layer conductor-backed structure (different thicknesses).
% Layer 1: er = 4 - 0.01j, ur = 1, thk = 1 mm
er2 = 4 - 0.1j;
ur2 = [];
thk2 = 2.^linspace(-3, -0, 7);

nLayer.printStructure(er2, ur2, thk2(1), Title="Case 2");
tic;
gam2 = zeros(length(f), length(thk2));
for ii = 1:length(thk2)
    gam2(:, ii) = NL.calculate(f, er2, ur2, thk2(ii));
end
fprintf("Various thicknesses: ");
toc;

% Plot
figure;
plot(gam2, ".-", Linewidth=1);
hold on;
title("Case 2");
zplane([]);
legend(compose("t = %.2f mm", abs(thk2)));
grid on;

%% Calculate structure 3
% One layer conductor-backed resonant structure (different losses).
% Layer 1: er = 4 - 0.0001j, ur = 1, thk = 2 mm
er3 = 4 - 1j.*10.^linspace(-1, -4, 7).';
ur3 = [];
thk3 = [2];

nLayer.printStructure(er3(1, :), ur3, thk3, Title="Case 3");
tic;
gam3 = zeros(length(fRes), length(er3));
for ii = 1:length(er3)
    gam3(:, ii) = NL.calculate(fRes, er3(ii), ur3, thk3);
end
fprintf("Resonant case various loss factors: ");
toc;

% Plot
figure;
plot(gam3, ".-", Linewidth=1);
hold on;
title("Case 3");
zplane([]);
legend(compose("\\epsilon'' = %.1e", abs(imag(er3))));
grid on;

%% Calculate structure 4
% One layer conductor-backed resonant structure (different conductivities).
% Layer 1: er = 4 - 0j, ur = 1, thk = 2 mm, sigma = 1e8
er4 = 4;
ur4 = [];
thk4 = [2];
sigma4 = 10.^linspace(1, 5, 9).';

nLayer.printStructure(er4(1, :), ur4, thk4, BackingConductivity=sigma4(1), ...
    Title="Case 4");
tic;
gam4 = zeros(length(fRes), length(sigma4));
for ii = 1:length(sigma4)
    gam4(:, ii) = NL.calculate(fRes, er4, ur4, thk4, ...
        BackingConductivity=sigma4(ii));
end
fprintf("Resonant case various conductivities: ");
toc;

% Plot
figure;
plot(gam4, ".-", Linewidth=1);
hold on;
title("Case 4");
zplane([]);
legend(compose("\\sigma = %.1e S/mm", abs(sigma4)));
grid on;
