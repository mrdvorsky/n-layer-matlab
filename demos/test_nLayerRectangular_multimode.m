clc;
clear;
close all;

%% Inputs
clc;
clear;
close all;

%% Inputs
wgR = 15.8;
f = linspace(20, 30, 51).';

numModes = 10;
numModes_show = 2;

%% Create nLayer Object
NL = nLayerRectangular(5, 4, Band="x");

%% Calculate structure 1
% One layer conductor-backed structure (different losses).
% Layer 1: er = 4 - 0.0001j, ur = 1, thk = 2 mm
er1 = [2 - 0.3j];
ur1 = [];
thk1 = [10];

NL.printStructure(er1, ur1, thk1, Title="Case 1");
gam1 = NL.calculate(f, er1, ur1, thk1);
gam1_modes = gam1(:, 1:numModes_show, 1:numModes_show);

legend_i = repmat(1:1:numModes_show, numModes_show, 1);
legend_j = legend_i.';
legendLabels = compose("S%d%d", legend_i(:), legend_j(:));

% Plot
figure;
plot(gam1_modes(:, :), ".-", Linewidth=1);
% plot(gam1, ".-", Linewidth=1);
hold on;
title("Case 1");
zplane([]);
legend(legendLabels);
grid on;


