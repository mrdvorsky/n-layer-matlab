clc;
clear;
close all;

%% Inputs
er = [4.0 - 0.4j];
ur = [1];
thk = [1.2];

f = linspace(32, 40, 2001);
% f = 36.61;

wgR = 30/128 * 25.4;
numModes = 3;

%% nLayerRectangularOld
tic;
NL1 = nLayerCircularTE_old(numModes, waveguideR=wgR, ...
    convergenceAbsTol=1e-5, verbosity=1);
toc;

%% nLayerRectangularFast
tic;
NL2 = nLayerCircularTE_fast(numModes, waveguideR=wgR, ...
    verbosity=1);
toc;
NL2.waveguideEr = 1;

%% Calculate
tic;
gam1 = NL1.calculate(f, er, ur, thk);
toc;

tic;
gam2 = NL2.calculate(f, er, ur, thk);
toc;

relErr = abs(max(gam1 - gam2))

%% Plot
figure;
plot(gam1, "-", LineWidth=1.5);
hold on;
plot(gam2, "-", LineWidth=1.5);
zplane([]);
grid on;
legend(["Original", "New"]);

% figure;
% plot(f, real(gam2), "-", LineWidth=1.5);
% hold on;
% plot(f, imag(gam2), "-", LineWidth=1.5);
% grid on;
% legend(["Real", "Imag"]);


% nLayer.plotModeStruct(NL2.modeStructs{1}, SizeX=3*wgR, SizeY=3*wgR)



