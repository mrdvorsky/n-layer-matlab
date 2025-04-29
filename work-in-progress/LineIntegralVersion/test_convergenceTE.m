clc;
clear;
close all;

%% Inputs
m = 3;
n = 3;

NL = nLayerRectangular(m, n, waveguideBand="Ka", modeSymmetryX="None", modeSymmetryY="None");
wgA = NL.waveguideA;
wgB = NL.waveguideB;

mode1 = NL.modeStructs(find(strcmp(NL.modeLabels, "TE_{1,2}")));
mode2 = NL.modeStructs(find(strcmp(NL.modeLabels, "TE_{1,2}")));

kc1 = mode1.CutoffWavenumber;
kc2 = mode2.CutoffWavenumber;

L = 1;
kphi(1, :) = linspace(0, 2*pi, 401);
[kr, kr_w] = fejer2_halfOpen(16000, L);

%% Calculate Integrals
kx = kr .* cos(kphi);
ky = kr .* sin(kphi);

spec1p = cos(kphi) .* mode1.ExSpec(kx, ky) + sin(kphi) .* mode1.EySpec(kx, ky);
spec2p = cos(kphi) .* mode2.ExSpec(kx, ky) + sin(kphi) .* mode2.EySpec(kx, ky);

spec1m = sin(kphi) .* mode1.ExSpec(kx, ky) - cos(kphi) .* mode1.EySpec(kx, ky);
spec2m = sin(kphi) .* mode2.ExSpec(kx, ky) - cos(kphi) .* mode2.EySpec(kx, ky);

intp = sum(spec1p .* spec2p .* kr .* kr_w ./ (L + kr), 1);
intm = sum(spec1m .* spec2m .* kr .* kr_w .* (L + kr), 1);


2 * trapz(kphi, intp + intm)

%% Plotting
figure;
plot(1:numel(kphi), intm, "", LineWidth=1.5);
grid on;

figure;
plot(1:numel(kphi), intp, "", LineWidth=1.5);
grid on;


