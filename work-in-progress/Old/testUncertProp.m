clc;
clear;
close all;

%% Inputs
f(:, 1) = linspace(32, 40, 41);

er = 25 - 0.01j;
thk(1, :) = linspace(0.01, 0.6, 61);

noiseStd = 10.^(-40 ./ 20);

%% nLayer
NL = nLayerCircularTE(3, waveguideR=5.8);
% NL = nLayerRectangular(3, 2, waveguideBand="ka");
NL.convergenceAbsTol = 1e-6;

NLsolver = nLayerInverse(1);
NLsolver.setLayersToSolve(Erp=[1], Erpp=[1], Thk=[1], Urp=[], Urpp=[]);
NLsolver.setInitialValues(Er=er);

erUncert = zeros(size(thk));
thkUncert = zeros(size(thk));
for tt = 1:numel(thk)
    tt
    NLsolver.setInitialValues(Thk=thk(tt));
    Uncert = NLsolver.computeParameterUncertainty(NL, f, NoiseStd=noiseStd);

    erUncert(tt) = Uncert.er;
    thkUncert(tt) = Uncert.thk;
end

%% Plot
figure;
semilogy(thk, thkUncert, "", LineWidth=1.5);
grid on;
xlabel("thk, mm");
ylabel("thk Uncertainty");
title("thk");

figure;
semilogy(thk, real(erUncert), "", LineWidth=1.5);
grid on;
xlabel("thk, mm");
ylabel("er^{'} Uncertainty");
title("er^{'}");

figure;
semilogy(thk, imag(erUncert), "", LineWidth=1.5);
grid on;
xlabel("thk, mm");
ylabel("er^{''} Uncertainty");
title("er^{''}");




