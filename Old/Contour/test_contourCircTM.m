clc;
clear;
close all;

%% Inputs
wgR = 1;

f = 130;

er = 1 - 0.01j;
thk = 10;

kr = linspace(0, 300, 10001);

L = 1;

%% nLayer Object
NL = nLayerCircularTM(1, waveguideR=wgR);

gam1 = NL.calculate(f, er, [], thk);

%% Integral
kc = NL.modeCutoffs;
k0 = 2*pi * f ./ NL.speedOfLight;

modeInt = @(kr) kr.^3 .* besselj(0, wgR .* kr).^2 ./ (kr.^2 - kc.^2).^2;
gamInt = @(kr) nLayerCircularTM.computeGamma0(kr, k0, er, 1, thk);
int = @(kr) modeInt(kr) .* gamInt(kr);

tic;
Imn = integral(int, 0, inf);
toc;

tic;
Imn2_1 = integral(@(kr) int(kr * (1 + 1j)), 0, L);
Imn2_2 = integral(@(kr) int(kr + 1j*L), L, inf);
toc;


tic;
Imn3 = integral(@(kr) int(complex(kr, kr ./ (1 + kr))) ...
    .* (1 + 1j ./ (1 + kr).^2), 0, inf);
toc;

err_db2 = db(Imn - (Imn2_1*(1 + 1j) + Imn2_2))
err_db3 = db(Imn - Imn3)

%% Reflection Coefficient
k0m = conj(sqrt(k0.^2 - kc.^2));
A = kc .* Imn .* k0m .* besselj(1, kc);
B = 0.5 .* besselj(1, kc) .* kc .* k0;

gam2 = (B - A) ./ (B + A);

%% Paths
valStraight = int(kr);

krCurved = complex(kr, L * kr ./ (kr + L));
valCurved = int(krCurved);

%% Plotting
figure;
plot(kr, real(valStraight), "", LineWidth=1.5);
hold on;
plot(kr, imag(valStraight), ":", LineWidth=1.5);
xlabel("kr");
grid on;

figure;
plot(kr, real(valCurved), "", LineWidth=1.5);
hold on;
plot(kr, imag(valCurved), ":", LineWidth=1.5);
xlabel("kr");
grid on;

figure;
semilogx(kr, db(valCurved), "", LineWidth=1.5);
hold on;
semilogx(kr, db(valStraight), ":", LineWidth=1.5);
xlabel("kr");
grid on;





