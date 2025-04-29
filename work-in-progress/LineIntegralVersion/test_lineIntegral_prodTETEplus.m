clc;
clear;
close all;

%% Inputs
m = 3;
n = 3;

NL = nLayerRectangular(m, n, waveguideBand="K", modeSymmetryX="None", modeSymmetryY="None");
wgA = NL.waveguideA;
wgB = NL.waveguideB;

mode1 = NL.modeStructs(find(strcmp(NL.modeLabels, "TE_{1,1}")));
mode2 = NL.modeStructs(find(strcmp(NL.modeLabels, "TE_{3,3}")));

kc1 = mode1.CutoffWavenumber;
kc2 = mode2.CutoffWavenumber;

krMax = 600;

[kr(:, 1), kr_w(:, 1)] = fejer2(400, 0, krMax);
[kphi(1, :), kphi_w(1, :)] = fejer2(100, 0, pi/2);
kphi_w = 4 * kphi_w;

kx = kr .* cos(kphi);
ky = kr .* sin(kphi);

%% Analytical
spec1 = cos(kphi) .* mode1.ExSpec(kx, ky) + sin(kphi) .* mode1.EySpec(kx, ky);
spec2 = cos(kphi) .* mode2.ExSpec(-kx, -ky) + sin(kphi) .* mode2.EySpec(-kx, -ky);

specR1 = kr .* sum(kphi_w .* spec1 .* spec2, 2) ./ (1 + kr);

%% Line Integrals
m1xy = (mode1.Hertz_x + 1j*mode1.Hertz_y);
m2xy = (mode2.Hertz_x + 1j*mode2.Hertz_y).';
r12 = abs(m1xy - m2xy);
phi12_m = angle(m1xy - m2xy);
phin12_m = mode1.Hertz_ang - mode2.Hertz_ang.';
phin12_p = mode1.Hertz_ang + mode2.Hertz_ang.';
Az12 = (mode1.Hertz_Hz .* mode1.Hertz_w) .* (mode2.Hertz_Hz .* mode2.Hertz_w).';

[~, sortInd] = sort(r12(:));

r12 = r12(sortInd).';
phi12_m = phi12_m(sortInd).';
phin12_m = phin12_m(sortInd).';
phin12_p = phin12_p(sortInd).';
Az12 = Az12(sortInd).';

% tic;
% specR2 = -kr .* sum(Az12 ...
%     .* (besselj(0, r12 .* kr) .* cos(phin12_m) ...
%     -   besselj(2, r12 .* kr) .* cos(phin12_p - 2*phi12_m)), ...
%     2) .* pi ./ (1 + kr);
% toc;

%% Integrals
I1 = sum(specR1 .* kr_w)
% I2 = sum(specR2 .* kr_w)

tic;
[gInt0, gInt2] = getWeightingFunctions(krMax, min(r12(r12 > 0)), max(r12), 4000, 8000);
toc;

Az12tmp = Az12;
% Az12tmp(r12 == 0) = 0;
I3 = -sum(Az12tmp ...
    .* (gInt0(r12) .* cos(phin12_m) ...
    -   gInt2(r12) .* cos(phin12_p - 2*phi12_m)), ...
    2)

I13 = I1 - I3

Az12tmp(r12 == 0) = 0;
I3_ext = -sum(Az12tmp ...
    .* (gInt0(r12) .* cos(phin12_m) ...
    -   gInt2(r12) .* cos(phin12_p - 2*phi12_m)), ...
    2)

%% Plot Spectrum
% figure;
% plot(kr, specR1, "", LineWidth=1.5);
% grid on;
% 
% figure;
% plot(kr, specR2, "", LineWidth=1.5);
% grid on;
% 
% figure;
% plot(kr, db((specR1 - specR2) ./ max(abs(specR1))), "", LineWidth=1.5);
% grid on;

%% Plot Integral
% figure;
% plot(1:numel(r12), r12, "", LineWidth=1.5);
% grid on;
% 
% figure;
% plot(1:numel(r12), Az12tmp .* gInt0(r12) .* cos(phin12_m), "", LineWidth=1.5);
% grid on;
% 
% figure;
% plot(1:numel(r12), Az12tmp .* gInt2(r12) .* cos(phin12_p - 2*phi12_m), "", LineWidth=1.5);
% grid on;

