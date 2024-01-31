clc;
clear;
close all;

%% Inputs
er = [4.0 - 0.4j];
ur = [1];
thk = [1.2];

f = linspace(26.5, 40, 5001);

wgR = 4.5;

%% nLayerCircularTM_old
tic;
NL1 = nLayerCircularTM_old(3, waveguideR=wgR, ...
    convergenceAbsTol=1e-4, verbosity=1);
toc;

%% nLayerRectangularFast
tic;
NL2 = nLayerCircularTM_fast(3, waveguideR=wgR, ...
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
% 
% figure;
% plot(f, real(gam2), "-", LineWidth=1.5);
% hold on;
% plot(f, imag(gam2), "-", LineWidth=1.5);
% grid on;
% legend(["Real", "Imag"]);


% nLayer.plotModeStruct(NL2.modeStructs{1}, SizeX=3*wgR, SizeY=3*wgR)



