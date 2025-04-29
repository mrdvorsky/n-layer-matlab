clc;
clear;
close all;

%% Inputs
Nm(1, 1) = 2;
Nrho(1, 1) = 1024;
L(1, 1) = 1;
Lc(1, 1) = 1;
Lch(1, 1) = 2;
Lcw(1, 1) = 10;

%% Helper
krToKrc = @(x) x + 1j * (Lch.*Lcw) * x ./ sqrt((Lch + x.^2).*(3*Lcw.^2 + x.^2));
krToKrcPrime = @(x) 1 + 1j * (Lch.*Lcw) * (3*Lch*Lcw.^2 - x.^4) ./ ((Lch + x.^2).*(3*Lcw.^2 + x.^2)).^1.5;

krcToKr = @(x) x - 1j * (Lch.*Lcw) * real(x) ./ sqrt((Lch + real(x).^2).*(3*Lcw.^2 + real(x).^2));
krcToKrPrime = @(x) 1 - 1j * (Lch.*Lcw) * (3*Lch*Lcw.^2 - real(x).^4) ./ ((Lch + real(x).^2).*(3*Lcw.^2 + real(x).^2)).^1.5;

%% Calculate Weights
[kr, kr_weights] = fejer2_halfOpen(Nrho, L);
krc = krToKrc(kr);

moment_weights = 1 ...
    .* sin(2*(Nm) .* acot(sqrt(kr./Lc))) ...
    ./ sin(2 * acot(sqrt(kr./Lc))) - Nm;

%% Plotting
figure;
plot(1:numel(kr), real(moment_weights), "", LineWidth=1.5);
hold on;
plot(1:numel(kr), imag(moment_weights), "", LineWidth=1.5);
grid on;






