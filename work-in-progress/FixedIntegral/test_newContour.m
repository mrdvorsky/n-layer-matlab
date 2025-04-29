clc;
clear;
% close all;

%% Inputs
f = @(x) 1 ./ (1 + x).^2;

Lc = 1;

%% Helper
krToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + (x));
krToKrcPrime = @(x) 1 + 1j .* Lc.^2 ./ ((x) + Lc).^2;

ang = deg2rad(270);
krcToKr1 = @(x) 0.5 * (exp(-0.5j*ang).*sqrt(exp(1j*ang).*(x.^2 + Lc*(2-2j)*x + 2j*Lc.^2)) + x - Lc - 1j*Lc);
krcToKr2 = @(x) 0.5 * (-exp(-0.5j*ang).*sqrt(exp(1j*ang).*(x.^2 + Lc*(2-2j)*x + 2j*Lc.^2)) + x - Lc - 1j*Lc);

% signReal = @(x) sign((krcError1(x) < krcError2(x)) - 0.5);
% krcToKr = @(x) 0.5 * (signReal(x) .* sqrt(x.^2 + Lc*(2-2j)*x + 2j*Lc.^2) + x - Lc - 1j*Lc);

signReal34 = @(x) real(x) < -imag(x);
% krcToKr3 = @(x) krcToKr1(x).*(signReal34(x)) + krcToKr2(x).*(~signReal34(x));
% krcToKr4 = @(x) krcToKr1(x).*(~signReal34(x)) + krcToKr2(x).*(signReal34(x));

% krcError3 = @(x) abs(x - krcToKr3(krToKrc(x)));
% krcError4 = @(x) abs(x - krcToKr4(krToKrc(x)));

% krcToKr = @(x) krcToKr3(x) + krcToKr4(x);
% krcToKr = @(x) x - + 1j .* Lc .* real(x) ./ (Lc + (real(x)));

krcToKr = @(x) krcToKr1(x).*(signReal34(x)) + krcToKr2(x).*(~signReal34(x));


krcToKrPrime = @(x) 0.5 * (1 + signReal34(x) .* (x + Lc - 1j*Lc) ...
    ./ sqrt(x.^2 + Lc*(2-2j)*x + 2j*Lc.^2));


%% Integral
% I1 = integral(f, 0, inf);
% I2 = integral(@(x) f(krToKrc(x)) .* krToKrcPrime(x), 0, inf);
% I3 = integral(@(x) f(krcToKr(x)) .* krcToKrPrime(x), 0, inf);

%% Plotting
rPlot(:, 1) = linspace(-10, 10, 1000);
iPlot(1, :) = linspace(-10, 10, 1000);
cPlot = rPlot + 1j*iPlot;

figure;
showImage(rPlot, iPlot, krcToKr1(cPlot), DisplayFormat="MagPhase");
axis normal;
colormap hsv;
clim([-180, 180]);

figure;
showImage(rPlot, iPlot, krcToKr2(cPlot), DisplayFormat="MagPhase");
axis normal;
colormap hsv;
clim([-180, 180]);

figure;
showImage(rPlot, iPlot, krcToKr(cPlot), DisplayFormat="MagPhase");
axis normal;
colormap hsv;
clim([-180, 180]);

% figure;
% showImage(rPlot, iPlot, krcToKr4(cPlot), DisplayFormat="MagPhase");
% axis normal;
% colormap hsv;
% clim([-180, 180]);











