clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(0, 1, waveguideR=0.5, modeSymmetryAxial="TM");
f = 300;

%% Modes
specE1 = NL.waveguideModes.WeSpec;

wgR = NL.waveguideR;
kc = NL.mode_kc0;

%% Integral
% JH = @(x) besselj(0, x, 1).*besselh(0, 1, x, 1).*exp(1j*real(x));
JH = @(x) besselj(0, x).*besselh(0, 1, x);

intFun1 = @(kr) -(2*pi) * specE1(0, 0, kr, 0).^2 .* kr;
intFun2 = @(kr) 2*besselj(0, wgR.*kr).^2.*(kr./(kr.^2 - kc.^2)).^2 .* kr;
intFun2 = @(kr) 2*(JH(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;
% intFun2 = @(kr) 2*(besselj(0, wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;
% intFun2 = @(kr) 2*(besselh(0, wgR.*kr) .* besselh(0, 1, wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;

int1 = integral(intFun2, 0, 10000)
% int2 = integral(intFun2, 0, 10 + 10j) + integral(intFun2, 10 + 10j, 20)
int3 = integral(intFun2, 0, 500000*(1 + 1j))
int4 = integral(intFun2, 0, 500000*(-1 + 1j))
% int3 = integral(intFun2, 0, -50*(1 + 1j))
% int3 = integral(intFun3, 0.01, 30*(1 + 1j))

% int3 + Correction


%% Residue
int3_res = integral(intFun2, 4, 4, Waypoints=[6 + 1j, 6 - 1j])


%% Plot
x(:, 1) = linspace(0, 20, 300);
y(1, :) = linspace(-1, 10, 300);


figure;
showImage(x, y, intFun2(x + 1j*y), DisplayFormat="Magnitude");
clim([0, 0.1])

%% Line Plot
% kr = linspace(0, 50, 1000) * (1 + 1j);
kr = linspace(-50, 50, 1000);

figure;
plot(real(kr), real(intFun1(kr)), "", LineWidth=1.5);
ylim(ylim());
hold on;
plot(real(kr), real(intFun2(kr)), ".", LineWidth=1.5);
grid on;

figure;
plot(real(kr), imag(intFun1(kr)), "", LineWidth=1.5);
ylim(ylim());
hold on;
plot(real(kr), imag(intFun2(kr)), ".", LineWidth=1.5);
grid on;






