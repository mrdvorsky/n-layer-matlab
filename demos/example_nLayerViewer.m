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
er1_start = [2.1];
ur1_start = [];
thk1_start = [3];

NL1 = nLayerRectangular(3, 2, waveguideBand="Ka");

figure;
nLayerViewer(er1_start, ur1_start, thk1_start, NL1, f1);

%% Example 2: Compare 6 modes to 1 mode
f2 = linspace(26.5, 40, 41).';
er2_start = [1, 2.1 - 0.001j];
ur2_start = [];
thk2_start = [0.01, 3];

NL2_1 = nLayerRectangular(3, 2, waveguideBand="Ka");
NL2_2 = nLayerRectangular(1, 0, waveguideBand="Ka");

figure;
nLayerViewer(er2_start, ur2_start, thk2_start, NL2_1, f2, NL2_2, f2, ...
    Legend=["6 modes", "1 mode"]);



