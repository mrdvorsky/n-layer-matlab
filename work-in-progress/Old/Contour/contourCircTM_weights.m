clc;
clear;
close all;

%% Inputs
wgR = 1;
f = 130;

er = [4 - 0.1j, 10];
thk = [1, 10];

L = pi;
intOrder = 32;

Lc = 1.5;

%% nLayer Object
NL = nLayerCircularTM(1, waveguideR=wgR);

%% Setup Integral
kc = NL.modeCutoffs;
k0 = 2*pi * f ./ NL.speedOfLight;

modeInt = @(kr) kr.^3 .* besselj(0, wgR .* kr).^2 ./ (kr.^2 - kc.^2).^2;
gamInt = @(kr) nLayerCircularTM.computeGamma0(kr, k0, er, 1, thk);
int = @(kr) modeInt(kr) .* gamInt(kr);

%% Integration Paths
xToKr = @(x) L * (1 - x) ./ x;
xToKr_weights = @(x) L ./ (x.^2);

% Curved
xToKrc = @(x) xToKr(x) + 1j .* Lc .* xToKr(x) ./ (Lc + xToKr(x));
xToKrc_weights = @(x) xToKr_weights(x) .* (1 + 1j .* Lc.^2 ./ (xToKr(x) + Lc).^2);

%% General Integral
Imn1 = integral(int, 0, inf, RelTol=1e-8, AbsTol=0);

%% Fejer Quadrature
[x, x_weights] = clenshawCurtis(intOrder, 0.0001, 1);

Imn2 = sum(int(xToKr(x)) .* xToKr_weights(x) .* x_weights);

%% Fejer Quadrature Complex Path
Imn3 = sum(int(xToKrc(x)) .* xToKrc_weights(x) .* x_weights);

%% Weighted Quadrature
wFun4 = @(x) modeInt(xToKr(x)) .* xToKr_weights(x);
fFun4 = @(x) gamInt(xToKr(x));

[x4, x4_weights] = clenshawCurtis(intOrder, 0.0001, 1, WeightingFunction=wFun4);

Imn4 = sum(fFun4(x4) .* x4_weights);
% Imn4 = integral(@(x) wFun4(x) .* fFun4(x), 0.0001, 1, RelTol=1e-8, AbsTol=0);

%% Weighted Quadrature Complex
wFun5 = @(x) modeInt(xToKrc(x)) .* xToKrc_weights(x);
fFun5 = @(x) gamInt(xToKrc(x));

[x5, x5_weights] = clenshawCurtis(intOrder, 0.0001, 1, WeightingFunction=wFun5);

Imn5 = sum(fFun5(x5) .* x5_weights);
% Imn5 = integral(@(x) wFun5(x) .* fFun5(x), 0.0001, 1, RelTol=1e-8, AbsTol=0);

%% Weight Calculation
test6 = @(x) wFun5(0.4999 .* cos(x) + 0.5) .* sin(x) .* cos((100)*x);

w6_1 = integral(test6, 0, pi);

[x6, x6_weights, x6_err] = fejer2(32000, 0, pi);

w6_2 = sum(x6_weights .* test6(x6));
w6_err = sum(x6_err .* test6(x6))

%% Error
% err_db2 = db(Imn1 - Imn2)
% err_db3 = db(Imn1 - Imn3)
% err_db4 = db(Imn1 - Imn4)
% err_db5 = db(Imn1 - Imn5)

%% Plotting
xPlot = linspace(0.0001, 1, 1001);

figure;
plot(xPlot, real(int(xToKr(xPlot)) .* xToKr_weights(xPlot)), "", LineWidth=1.5);
hold on;
plot(xPlot, imag(int(xToKr(xPlot)) .* xToKr_weights(xPlot)), "", LineWidth=1.5);
grid on;
title("Unweighted Real");

figure;
plot(xPlot, real(fFun4(xPlot)), "", LineWidth=1.5);
hold on;
plot(xPlot, imag(fFun4(xPlot)), "", LineWidth=1.5);
grid on;
title("Weighted Real");

figure;
plot(xPlot, real(int(xToKrc(xPlot)) .* xToKrc_weights(xPlot)), "", LineWidth=1.5);
hold on;
plot(xPlot, imag(int(xToKrc(xPlot)) .* xToKrc_weights(xPlot)), "", LineWidth=1.5);
grid on;
title("Unweighted Complex");

figure;
plot(xPlot, real(fFun5(xPlot)), "", LineWidth=1.5);
hold on;
plot(xPlot, imag(fFun5(xPlot)), "", LineWidth=1.5);
grid on;
title("Weighted Complex");

%% Plot Weights
figure;
plot(x5, real(wFun5(x5)), "", LineWidth=1.5);
hold on;
plot(x5, imag(wFun5(x5)), "", LineWidth=1.5);

figure;
plot(x5, real(x5_weights), "", LineWidth=1.5);
hold on;
plot(x5, imag(x5_weights), "", LineWidth=1.5);



