clc;
clear;
close all;

%% Inputs
wgR = 1;

x(:, 1) = 50.1 * linspace(-1, 1, 2001);
y(1, :) = 50.1 * linspace(-1, 1, 2001);

m = 1;
n = 1;

kc = besselj_der_zeros(n, m) ./ wgR;

thetaPlot(1, :) = [0, 90];
rhoPlot(:, 1) = linspace(0, 20, 101);

%% Coordinates
r = hypot(x, y);
phi = angle(x + 1j*y);

%% Calculate Fields
modeWin = (r <= wgR);

Er = n .* besselj(n, kc * r) .* sin(n .* phi) ./ r .* modeWin;
Ephi = kc .* besselj_der(n, kc * r) .* cos(n .* phi) .* modeWin;
Er(r == 0) = 0.5*(n == 1);

Ex = Er .* cos(phi) - Ephi .* sin(phi);
Ey = Er .* sin(phi) + Ephi .* cos(phi);

%% Spectrums
[kx, ky] = fftCoordinates(x, y, ApplyFftShift=true);
kr = hypot(kx, ky);
kphi = angle(kx + 1j*ky);

Ex_hat = real(fftshift(fft2(ifftshift(Ex))));
Ey_hat = real(fftshift(fft2(ifftshift(Ey))));

%% Interp
Ex_hat_int = griddedInterpolant({kx, ky}, Ex_hat);
Ey_hat_int = griddedInterpolant({kx, ky}, Ey_hat);

%% Plot
kxPlot = rhoPlot.*cos(thetaPlot);
kyPlot = rhoPlot.*sin(thetaPlot);

figure;
plot(rhoPlot, Ex_hat_int(kxPlot, kyPlot), "", LineWidth=1.5);
% plot(rhoPlot, Ey_hat_int(kxPlot, kyPlot), "", LineWidth=1.5);
grid on;

figure;
% plot(rhoPlot, Ex_hat_int(kxPlot, kyPlot), "", LineWidth=1.5);
plot(rhoPlot, Ey_hat_int(kxPlot, kyPlot), "", LineWidth=1.5);
grid on;

fit_Ey = Ey_hat_int(kxPlot, kyPlot);
fit_Ey = fit_Ey(:, 1);
fit_k = rhoPlot;


%% Plotting
scale = max(abs(cat(3, Ex_hat, Ey_hat)), [], "all");

figure;
subplot(1, 2, 1);
showImage(kx, ky, Ex_hat, DisplayFormat="Magnitude");
clim(scale * [0, 1]);
title("E_{x}");

subplot(1, 2, 2);
showImage(kx, ky, Ey_hat, DisplayFormat="Magnitude");
clim(scale * [0, 1]);
title("E_{y}");


%% Bessel Function Zeros
function [x, x_init] = besselj_zeros(v, m)
    x_init = pi * (m + 0.5*v - 0.25);
    x = fzero(@(y) besselj(v, y), x_init);
end

function [y] = besselj_der(v, x)
    y = 0.5 * (besselj(v - 1, x) - besselj(v + 1, x));
end

function [x, x_init] = besselj_der_zeros(v, m)
    x_init = pi * (m + 0.5*v + 0.25 - (v ~= 0));
    x = fzero(@(y) besselj_der(v, y), x_init);
end

function [y] = besselInt(kr, n, kc, wgR)
    y = wgR .* (kc .* besselj_der(n, wgR .* kr) .* besselj_der(n - 1, wgR .* kc) ...
        - kr .* besselj_der(n - 1, wgR .* kr) .* besselj_der(n, wgR .* kc)) ...
        ./ (kr.^2 - kc.^2);
end

