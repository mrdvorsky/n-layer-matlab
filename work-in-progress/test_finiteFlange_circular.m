clc;
clear;
close all;

%% Inputs
finiteFlange_r = 15;
finiteFlange_ro1 = 18;
finiteFlange_ro2 = 19;

numModes = 10;
numModes_flange = 10;

f(:, 1) = linspace(32, 40, 400);
er = {4 - 0.004j};
ur = {1};
thk = {2};

%% nLayer
tic;
NLcirc = nLayerCircular(0, numModes, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
NLcirc.calculate(f(1), er, ur, thk);
toc;

NLflange1 = nLayerCoaxial(0, numModes_flange, modeSymmetryAxial="TE", ...
    waveguideRi=finiteFlange_r, waveguideRo=finiteFlange_ro1);
NLflange2 = nLayerCoaxial(0, numModes_flange, modeSymmetryAxial="TE", ...
    waveguideRi=finiteFlange_r, waveguideRo=finiteFlange_ro1);

tic;
NL1 = nLayerOpenEnded([NLcirc.waveguideModes, NLflange1.waveguideModes], ...
    frequencyRange=f);
NL2 = nLayerOpenEnded([NLcirc.waveguideModes, NLflange2.waveguideModes], ...
    frequencyRange=f);
NL1.calculate(f(1), er, ur, thk);
NL2.calculate(f(1), er, ur, thk);
toc;

%% Plot
gamCirc = NLcirc.calculate(f, er, ur, thk);
gam1 = NL1.calculate(f, er, ur, thk);
gam2 = NL2.calculate(f, er, ur, thk);

figure;
plot(f, db([gamCirc, gam1, gam2]), "", LineWidth=1.5);
grid on;

figure;
plot(f, rad2deg(angle([gamCirc, gam1, gam2])), "", LineWidth=1.5);
grid on;




% figure;
% nLayerViewer(er, ur, thk, NLcirc, f, NL1, f, NL2, f, ...
%     NumFrequencySamplesPerMarker=40);





