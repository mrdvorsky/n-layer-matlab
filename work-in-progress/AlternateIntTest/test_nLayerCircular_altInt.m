clc;
clear;
close all;

%% Inputs
f = 35;

er = {4 - 0.2j};
ur = {1};
thk = {1.2};

numModes = 1;

%% nLayer
NL = nLayerCircular(0, numModes, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
NL2 = nLayerCircularTE_old(numModes, waveguideR=NL.waveguideR, convergenceAbsTol=1e-4);
k0 = 2*pi .* f ./ NL.speedOfLight;

%% Calculate
W = @(kr) NL.modeStructs(1).EySpec(0, 0, kr, 0).^2 .* kr .* (2*pi);

beta = sqrt(k0.^2 - NL.modeStructs(1).CutoffWavenumber.^2);
Ke = sqrt(k0 ./ beta);

shouldBeOne = integral(W, 0, inf)
Ih = integral(@(kr) W(kr) .* gam0H(kr, k0, er, ur, thk), 0, inf);

gam1 = (1 - Ke.^2 .* Ih) ./ (1 + Ke.^2 .* Ih)



% Ih2 = integral(@(kr) W(kr) .* gam0H(kr, k0, er, ur, thk))


%% Old
gam2 = NL2.calculate(f, er, ur, thk)

%% Plotting
krPlot = linspace(0, 10, 1000);

figure;
plot(krPlot, real(gam1H(krPlot, k0, er, ur, thk)), "", LineWidth=1.5);
hold on;
plot(krPlot, imag(gam1H(krPlot, k0, er, ur, thk)), "", LineWidth=1.5);

figure;
plot(krPlot, abs(gam1H(krPlot, k0, er, ur, thk)), "", LineWidth=1.5);










