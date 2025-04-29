clc;
clear;
close all;

%% Inputs
m = 3;
n = 3;

NL = nLayerRectangular(m, n, waveguideBand="Ka", modeSymmetryX="None", modeSymmetryY="None");
wgA = NL.waveguideA;
wgB = NL.waveguideB;

mode1 = NL.modeStructs(find(strcmp(NL.modeLabels, "TM_{1,1}")));
mode2 = NL.modeStructs(find(strcmp(NL.modeLabels, "TM_{1,1}")));

kc1 = mode1.CutoffWavenumber;
kc2 = mode2.CutoffWavenumber;

L = 1;
kphi(1, :) = deg2rad(45);
[kr, kr_w] = fejer2_halfOpen(16000, L);

%% Analytical
kx = kr .* cos(kphi);
ky = kr .* sin(kphi);

spec1 = cos(kphi) .* mode1.ExSpec(kx, ky) + sin(kphi) .* mode1.EySpec(kx, ky);
spec2 = cos(kphi) .* mode2.ExSpec(kx, ky) + sin(kphi) .* mode2.EySpec(kx, ky);

%% Spectrums Line
Ez1L = mode1.Line_Ez;

d1L = mode1.Line_d;
x1L = mode1.Line_x;
y1L = mode1.Line_y;
ang1L = mode1.Line_ang;

[~, s1L, c1L] = polyFourierTransform2(-0.5*d1L, 0.5*d1L, Ez1L);

tic;
spec1_line = 0;
for nn = 1:size(c1L, 2)
    spec1_line = spec1_line + sum(c1L(:, nn).' ...
        .* exp(-1j .* (cos(kphi).*x1L + sin(kphi).*y1L).' .* kr) ...
        .* besselj_spherical(nn - 1, (s1L(:) .* cos(ang1L - kphi)).' .* kr), 2);
end
toc;

spec1_line = spec1_line .* kr ./ (kr.^2 - kc1.^2);



Ez2L = mode2.Line_Ez;

d2L = mode2.Line_d;
x2L = mode2.Line_x;
y2L = mode2.Line_y;
ang2L = mode2.Line_ang;

[~, s2L, c2L] = polyFourierTransform2(-0.5*d2L, 0.5*d2L, Ez2L);

tic;
spec2_line = 0;
for nn = 1:size(c2L, 2)
    spec2_line = spec2_line + sum(c2L(:, nn).' ...
        .* exp(-1j .* (cos(kphi).*x2L + sin(kphi).*y2L).' .* kr) ...
        .* besselj_spherical(nn - 1, (s2L(:) .* cos(ang2L - kphi)).' .* kr), 2);
end
toc;

spec2_line = spec2_line .* kr ./ (kr.^2 - kc2.^2);

%% Spec12_line
re1 = (cos(kphi).*x1L + sin(kphi).*y1L).';
re2 = (cos(kphi).*x2L + sin(kphi).*y2L).';

rj1 = (s1L(:) .* cos(ang1L - kphi)).';
rj2 = (s2L(:) .* cos(ang2L - kphi)).';

spec12_line = 0;
for nn = 1:size(c1L, 2)
    for mm = 1:size(c2L, 2)
        spec12_line = spec12_line + sum(c2L(:, nn).' ...
            .* exp(-1j .* (cos(kphi).*x2L + sin(kphi).*y2L).' .* kr) ...
            .* besselj_spherical(nn - 1, (s2L(:) .* cos(ang2L - kphi)).' .* kr), 2);
    end
end

%% Integrals
I1 = sum(spec1 .* spec2 .* kr .* kr_w ./ (L + kr))

I2 = sum(spec1_line .* spec2_line .* kr .* kr_w ./ (L + kr))

%% Plot Spectrum
figure;
plot(1:numel(kr), real(spec1), "", LineWidth=1.5);
hold on;
plot(1:numel(kr), imag(spec1), "", LineWidth=1.5);
grid on;

figure;
plot(1:numel(kr), real(spec1_line), "", LineWidth=1.5);
hold on;
plot(1:numel(kr), imag(spec1_line), "", LineWidth=1.5);
grid on;


figure;
plot(1:numel(kr), db((spec1_line - spec1) ./ max(abs(spec1))), "", LineWidth=1.5);
grid on;
ylim([-300, 0]);


