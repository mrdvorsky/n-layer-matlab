clc;
clear;
close all;

%% Inputs
m = 3;
n = 3;

NL = nLayerRectangular(m, n, waveguideBand="Ka", modeSymmetryX="None", modeSymmetryY="None");
wgA = NL.waveguideA;
wgB = NL.waveguideB;

mode1 = NL.modeStructs(find(strcmp(NL.modeLabels, "TE_{1,1}")));
mode2 = NL.modeStructs(find(strcmp(NL.modeLabels, "TE_{1,3}")));

kc1 = mode1.CutoffWavenumber;
kc2 = mode2.CutoffWavenumber;

[kr(:, 1), kr_w(:, 1)] = fejer2(200, 0, 10);
[kphi(1, :), kphi_w(1, :)] = fejer2(100, 0, pi/2);
kphi_w = 4 * kphi_w;

kx = kr .* cos(kphi);
ky = kr .* sin(kphi);

%% Analytical
spec1 = sin(kphi) .* mode1.ExSpec(kx, ky) - cos(kphi) .* mode1.EySpec(kx, ky);
spec2 = sin(kphi) .* mode2.ExSpec(kx, ky) - cos(kphi) .* mode2.EySpec(kx, ky);

specR1 = kr .* sum(kphi_w .* spec1 .* spec2, 2);

%% Line Integrals
m1xy = (mode1.Hertz_x + 1j*mode1.Hertz_y);
m2xy = (mode2.Hertz_x + 1j*mode2.Hertz_y).';
r12 = abs(m1xy - m2xy);
phi12_m = angle(m1xy - m2xy);
phin12_m = mode1.Hertz_ang - mode2.Hertz_ang.';
phin12_p = mode1.Hertz_ang + mode2.Hertz_ang.';
Az12 = (mode1.Hertz_Hz .* mode1.Hertz_w) .* (mode2.Hertz_Hz .* mode2.Hertz_w).';

r12 = r12(:).';
phi12_m = phi12_m(:).';
phin12_m = phin12_m(:).';
phin12_p = phin12_p(:).';
Az12 = Az12(:).';

tic;
specR2 = kr .* sum(Az12 ...
    .* (besselj(0, r12 .* kr) .* cos(phin12_m) ...
    +   besselj(2, r12 .* kr) .* cos(phin12_p - 2*phi12_m)), ...
    2) .* pi .* (kc1 .* kc2).^2 ./ ((kr.^2 - kc1.^2) .* (kr.^2 - kc2.^2));
toc;


% (kc1 .* kc2).^2 ./ ((kr.^2 - kc1.^2) .* (kr.^2 - kc2.^2));

%% Plot Spectrum
figure;
plot(kr, specR1, "", LineWidth=1.5);
grid on;

figure;
plot(kr, specR2, "", LineWidth=1.5);
grid on;

figure;
plot(kr, db((specR1 - specR2) ./ max(abs(specR1))), "", LineWidth=1.5);
grid on;




