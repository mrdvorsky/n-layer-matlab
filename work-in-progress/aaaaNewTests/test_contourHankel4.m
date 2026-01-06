clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(0, 1, waveguideR=0.5, modeSymmetryAxial="TM");
f = 300;

%% Modes
wgR = NL.waveguideR;
kc = NL.mode_kc0;

%% Integrands
JH1 = @(x) besselj(0, x).^2;
% JH2 = @(x) besselj(0, x, 1).*besselh(0, 1, x, 1).*exp(1j*real(x));
JH2 = @(x) besselj(0, x).*besselh(0, 1, x);

intFun1 = @(kr) 10*(JH1(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;
intFun2 = @(kr) 10*(JH2(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;

%% Full Integral
sqrt_bc = @(x) real(sqrt(x)) - 1j*abs(imag(sqrt(x)));
g = @(kr) 1 ./ sqrt_bc((2*pi - 0.001j).^2 - kr.^2);

int1 = integral(@(kr) g(kr) .* intFun1(kr), 0, 5000)
int2 = integral(@(kr) g(kr) .* intFun2(kr), 0, 5000*(1 + 1j))



int_res2 = integral(@(kr) g(kr) .* intFun2(kr), 4, 4, Waypoints=[5 + 1j, 5 - 1j])
% int_res3 = integral(@(kr) g(kr) .* intFun2(kr), 4, 4, Waypoints=[5 + 1j, 5 - 1j])





