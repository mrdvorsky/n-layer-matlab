clc;
clear;
close all;

%% Inputs
f = 35;

er = {4 - 0.8j};
ur = {1};
thk = {0.8};

numModes = 1;

%% nLayer
NL = nLayerCircular(0, numModes, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
NL2 = nLayerCircularTE_old(numModes, waveguideR=NL.waveguideR, convergenceAbsTol=1e-4);
k0 = 2*pi .* f ./ NL.speedOfLight;

beta = sqrt(k0.^2 - NL.modeStructs(1).CutoffWavenumber.^2);
Ke = sqrt(k0 ./ beta);

%% Calculate
W = @(kr) NL.modeStructs(1).EySpec(0, 0, kr, 0).^2 .* kr .* (2*pi);

shouldBeOne = integral(W, 0, inf)

intFun1 = @(kr) W(kr) - Ke.^2 .* W(kr) .* gam0H(kr, k0, er, ur, thk);
intFun2 = @(kr) W(kr) + Ke.^2 .* W(kr) .* gam0H(kr, k0, er, ur, thk);
Ih_1 = integral(intFun1, 0, inf);
Ih_2 = integral(intFun2, 0, inf);
gam1 = Ih_1 ./ Ih_2

%% Old
gam2 = NL2.calculate(f, er, ur, thk)

%% Plotting
krp = linspace(0, 2, 1000);

figure;
plot(krp, real(intFun1(krp)), "", LineWidth=1.5);
hold on;
plot(krp, imag(intFun1(krp)), "", LineWidth=1.5);
grid on;

hold on;
plot(krp, real(gam1 .* intFun2(krp)), ":", LineWidth=1.5);
hold on;
plot(krp, imag(gam1 .* intFun2(krp)), ":", LineWidth=1.5);
grid on;









