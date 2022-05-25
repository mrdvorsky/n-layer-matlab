% This example file shows how to use nLayerRectangular to
% calculate S11 of a rectangular waveguide looking into a multilayer
% structure.
%
% The basic usage of nLayerRectangular is as follows (example for x-band):
%   NL = nLayerRectangular(maxM, maxN, band="x");
%   S11 = NL.calculate(f, er, ur, thk);
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
wgBand = "x";                       % Use x-band waveguide dimensions
f = linspace(8.2, 12.4, 21).';      % Frequencies to calculate (column vector)
fRes = linspace(8.2, 12.4, 401).';  % Use more points for the resonant case

maxM = 3;                           % Max mode index for broad dimension
maxN = 2;                           % Max mode index for narrow dimension

convergenceTol = 1e-3;              % Tolerance for convergence (-60 dB)

%% Create nLayer Object
% Creating the object initially takes a short amount of time (~0.2 seconds)
% to precompute integral stuff. You can either create this once at the
% beginning of your program, or save the nLayer object to a file and load
% it (second line here). Waveguide dimensions can also be specified
% directly or using the waveguide band identifier. See documentation for
% nLayerRectangular for more details.

NL = nLayerRectangular(maxM, maxN, Band=wgBand);        % Create nLayer object.
% NL = load("nLayer_xBand").NL;                           % Optionally load nLayer object from file.
% NL = nLayerRectangular(maxM, maxN, A=22.86, B=10.16);   % Alternative way to specify waveguide dimensions.

NL.convergenceAbsTol = convergenceTol;  % Optionally specify tolerance for convergence.
NL.verbosity = 1;                       % Leave at 0 for no command line output.

%% Calculate structure 1
% Two layer conductor-backed structure (very low loss).
% Layer 1: er = 1 - 0.0001j, ur = 1, thk = 1 mm
% Layer 2: er = 2 - 0.0001j, ur = 1, thk = 5 mm
er1 = [1 - 0.0001j, 2 - 0.0001j];   % Row vector
ur1 = [];                           % Empty vector will be defaulted to all 1's
thk1 = [1, 5];                      % Row vector

NL.printStructure(er1, ur1, thk1, Title="Case 1");
tic;
% Call this function to calculate the reflection coefficient. The
% reflection coefficient will be calculated within the tolerance specified
% earlier. Optionally, a custom tolerance can be specified (second line).
gam1 = NL.calculate(f, er1, ur1, thk1);
fprintf("Two-layer low loss conductor backed: ");
toc;

% Plot
figure;
plot(gam1, "-", "Linewidth", 1);
hold on;
plot(gam1, ".", "Linewidth", 1.5);
title("Case 1");
zplane([]);
grid on;

%% Calculate structure 2
% Two layer infinte half-space structure
% Layer 1: er = 1 - 0.0001j, ur = 1, thk = 1 mm
% Layer 2: er = 2 - 0.0001j, ur = 1, thk = inf
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
plot(gam2, "-", "Linewidth", 1);
hold on;
plot(gam2, ".", "Linewidth", 1.5);
title("Case 2");
zplane([]);
grid on;

%% Calculate structure 3
% One layer infinte half-space structure with permeability
% Layer 1: er = 4 - 1j, ur = 3 - 0.1j, thk = inf
er3 = [4 - 1j];
ur3 = [3 - 0.1j];
thk3 = [inf];

NL.printStructure(er3, ur3, thk3, Title="Case 3");
tic;
gam3 = NL.calculate(f, er3, ur3, thk3);
fprintf("One-layer lossy half-space: ");
toc;

% Plot
figure;
plot(gam3, "-", "Linewidth", 1);
hold on;
plot(gam3, ".", "Linewidth", 1.5);
title("Case 3");
zplane([]);
grid on;

%% Calculate structure 4
% One-layer conductor backed resonant case
% Layer 1: er = 2 - 0.001j, ur = 1, thk = 10 mm
er4 = [2 - 0.001j];
ur4 = [];
thk4 = [10];

NL.printStructure(er4, ur4, thk4, Title="Case 4");
tic;
gam4 = NL.calculate(fRes, er4, ur4, thk4);
fprintf("One-layer conductor backed resonant: ");
toc;

% Plot
figure;
plot(gam4, "-", "Linewidth", 1);
hold on;
plot(gam4, ".", "Linewidth", 1.5);
title("Case 4");
zplane([]);
grid on;


%% Calculate structure 5
% Two-layer conductor backed frequency-variable structure
% Layer 1: er = 1 - 0.0001j to 1 - 1j, ur = 1, thk = 1 mm
% Layer 2: er = 2 - 0.0001j to 4 - 4j, ur = 1, thk = 5 mm
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
plot(gam5, "-", "Linewidth", 1);
hold on;
plot(gam5, ".", "Linewidth", 1.5);
title("Case 5");
zplane([]);
grid on;


