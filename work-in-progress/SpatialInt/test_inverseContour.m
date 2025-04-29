clc;
clear;
close all;

%% Inputs
x(:, 1) = linspace(-1, 30, 1000);
y(1, :) = linspace(-2, 2, 1000);

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

% krToKrc = @(x) x .* exp(0.25j*pi ./ (1 + x.^2));

%% Eval
alpha = krToKrc(x + 1j*y);

%% Plot
figure;
showImage(x, y, -alpha, DisplayFormat="Phase");
colormap jet;
hold on;
plot(x, imag(krToKrc(x)), "", LineWidth=1.5);
clim([-180, 180]);
grid on;
