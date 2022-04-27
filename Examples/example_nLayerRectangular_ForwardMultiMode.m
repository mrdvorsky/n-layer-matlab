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
f = linspace(20, 30, 201).';      % Frequencies to calculate (column vector)

maxM = 3;                           % Max mode index for broad dimension
maxN = 0;                           % Max mode index for narrow dimension

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
% Layer 1: er = 2 - 0.3j, ur = 1, thk = 10 mm
er1 = [2 - 0.3j];
ur1 = [];
thk1 = [10];

tic;
% Call this function to calculate the reflection coefficient. The
% reflection coefficient will be calculated within the tolerance specified
% earlier. Optionally, a custom tolerance can be specified (second line).
gam1 = NL.calculate(f, er1, ur1, thk1);
fprintf("Two-layer low loss conductor backed: ");
toc;

% Plot
figure;
plot(gam1(:, :), "-", "Linewidth", 1);
hold on;
plot(gam1(:, :), ".", "Linewidth", 1.5);
title("Case 1");
zplane([]);
grid on;

