clc;
clear;
close all;

%% Inputs
Nphi = 200;
Nkr = 8000;


Nkr2 = 4000;
Nr2 = 1000;

L = 10;
Mk = 20;

NL = nLayerRectangular(3, 3, waveguideBand="Ka", modeSymmetryX="None", modeSymmetryY="None");
mode1 = NL.modeStructs(find(strcmp(NL.modeLabels, "TM_{1,1}")));
mode2 = NL.modeStructs(find(strcmp(NL.modeLabels, "TM_{1,1}")));

%% Calculate
kc1 = mode1.CutoffWavenumber;
kc2 = mode2.CutoffWavenumber;

wgA = NL.waveguideA;
wgB = NL.waveguideB;

[kr(:, 1), kr_w(:, 1)] = fejer2_halfOpen(Nkr, L);
% [kr(:, 1), kr_w(:, 1)] = fejer2(Nkr, 0, 40);
[kphi(1, :), kphi_w(1, :)] = fejer2(Nphi, 0, pi/2);
kphi_w = 4 * kphi_w;


krToKrc = @(kr) kr + 1j * kr ./ (1 + kr.^2);
krToKrcPrime = @(kr) (1 + 1j .* (1 - kr.^2) ./ (1 + kr.^2).^2);

kx = krToKrc(kr) .* cos(kphi);
ky = krToKrc(kr) .* sin(kphi);

krc_w = kr_w .* krToKrcPrime(kr);

%% Analytical
spec1 = cos(kphi) .* mode1.ExSpec(kx, ky) + sin(kphi) .* mode1.EySpec(kx, ky);
spec2 = cos(kphi) .* mode2.ExSpec(kx, ky) + sin(kphi) .* mode2.EySpec(kx, ky);

specR1 = krToKrc(kr) .* sum(kphi_w .* spec1 .* spec2, 2) ./ (L + krToKrc(kr));

%% Line Integrals
r12 = hypot(mode1.Hertz_x - mode1.Hertz_x.', mode2.Hertz_y - mode2.Hertz_y.');
Az12 = (mode1.Hertz_Ez .* mode1.Hertz_w) .* (mode2.Hertz_Ez .* mode2.Hertz_w).';
r12 = r12(:).';
Az12 = Az12(:).';

tic;
rMax = max(r12);
r2(1, :) = linspace(0, rMax, Nr2);
[kr2, kr2_w] = fejer2_halfOpen(Nkr2, L);
% [kr2, kr2_w] = fejer2(Nkr2, 0, 1*3.14*L);
krc2 = krToKrc(kr2);

% kr2_w = kr2_w .* exp(-(kr2 ./ (13*L)).^2);

mom_w = sin(2*Mk * acot(sqrt((kr2)./L))) ...
     ./ sin(2    * acot(sqrt((kr2)./L)));

intVals = sum(mom_w .* besselj(0, r2 .* krToKrc(kr2)) .* kr2_w .* krToKrcPrime(kr2) ...
    .* krToKrc(kr2).^3 ./ ((krToKrc(kr2).^2 - kc1.^2) .* (krToKrc(kr2).^2 - kc2.^2)) .* (2*pi) ./ (L + krToKrc(kr2)), 1);
% for rr = 1:numel(r2)
%     rr
%     intVals(rr) = integral(@(kr) besselj(0, r2(rr) .* krToKrc(kr)) .* krToKrcPrime(kr) ...
%         .* krToKrc(kr).^3 ./ ((krToKrc(kr).^2 - kc1.^2) .* (krToKrc(kr).^2 - kc2.^2)) .* (2*pi) ./ (L + krToKrc(kr)), ...
%         0, inf, RelTol=1e-6);
% end

gr = griddedInterpolant(r2, intVals, "spline");

% tic;
% specR2 = sum(Az12 .* besselj(0, kr .* r12), 2) ...
%     .* kr.^3 ./ ((kr.^2 - kc1.^2) .* (kr.^2 - kc2.^2)) .* (2*pi) ./ (L + kr);
% toc;

I2 = sum(gr(r12) .* Az12)
toc;

%% Integrals
I1 = sum(specR1 .* krc_w)

errdb = db(I1 - I2)

%% Plot Spectrum
% figure;
% plot(1:numel(kr), specR1, "", LineWidth=1.5);
% grid on;
% 
% figure;
% plot(1:numel(kr), specR2, "", LineWidth=1.5);
% grid on;
% 
% figure;
% plot(1:numel(kr), db((specR1 - specR2) ./ max(abs(specR1))), "", LineWidth=1.5);
% grid on;




