clc;
clear;
close all;

%% Inputs
wgR = 1;

x(:, 1) = 50.1 * linspace(-1, 1, 2001);
y(1, :) = 50.1 * linspace(-1, 1, 2001);

m = 2;
n = 2;

kc = besselj_zeros(n, m) ./ wgR;

%% Coordinates
r = hypot(x, y);
phi = angle(x + 1j*y);

%% Calculate Fields
modeWin = (r <= wgR);

Er = -kc .* besselj_der(n, kc * r) .* cos(n .* phi) .* modeWin;
Ephi = n .* besselj(n, kc * r) .* sin(n .* phi) ./ r .* modeWin;
Ephi(r == 0) = 0.5*(n == 1);

Ex = Er .* cos(phi) - Ephi .* sin(phi);
Ey = Er .* sin(phi) + Ephi .* cos(phi);

%% Spectrums
[kx, ky] = fftCoordinates(x, y, ApplyFftShift=true);
kr = hypot(kx, ky);
kphi = angle(kx + 1j*ky);

scaleFactor = abs(x(1) - x(2)) * abs(y(1) - y(2));
Ex_hat = -1j*fftshift(fft2(ifftshift(Ex))) .* scaleFactor;
Ey_hat = -1j*fftshift(fft2(ifftshift(Ey))) .* scaleFactor;

Ex_hat2 = (2*pi) * kr .* cos(kphi) .* cos(n * kphi) .* besselInt(kr, n, kc, wgR);
Ey_hat2 = (2*pi) * kr .* sin(kphi) .* cos(n * kphi) .* besselInt(kr, n, kc, wgR);

errX_rel = abs(max(Ex_hat - Ex_hat2, [], "all") ./ max(Ex_hat2, [], "all"))
errY_rel = abs(max(Ey_hat - Ey_hat2, [], "all") ./ max(Ey_hat2, [], "all"))

%% Plotting
% scale = max(abs(cat(3, Ex, Ey, Er, Ephi)), [], "all");
% 
% figure;
% subplot(2, 2, 1);
% showImage(x, y, Er);
% clim(scale * [-1, 1]);
% title("E_{r}");
% 
% subplot(2, 2, 2);
% showImage(x, y, Ephi);
% clim(scale * [-1, 1]);
% title("E_{\phi}");
% 
% subplot(2, 2, 3);
% showImage(x, y, Ex);
% clim(scale * [-1, 1]);
% title("E_{x}");
% 
% subplot(2, 2, 4);
% showImage(x, y, Ey);
% clim(scale * [-1, 1]);
% title("E_{y}");

%% Plotting
scale = max(abs(cat(3, Ex_hat, Ey_hat, Ex_hat2, Ey_hat2)), [], "all");

figure;
subplot(2, 2, 1);
showImage(kx, ky, Ex_hat, DisplayFormat="Magnitude");
clim(scale * [0, 1]);
title("E_{x}");

subplot(2, 2, 2);
showImage(kx, ky, Ey_hat, DisplayFormat="Magnitude");
clim(scale * [0, 1]);
title("E_{y}");

scale2 = max(abs(cat(3, Ex_hat2, Ey_hat2)), [], "all");
subplot(2, 2, 3);
showImage(kx, ky, Ex_hat2, DisplayFormat="Magnitude");
clim(scale2 * [0, 1]);
title("E_{x}");

subplot(2, 2, 4);
showImage(kx, ky, Ey_hat2, DisplayFormat="Magnitude");
clim(scale2 * [0, 1]);
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
    y = wgR .* (kc .* besselj(n, wgR .* kr) .* besselj(n - 1, wgR .* kc) ...
        - kr .* besselj(n - 1, wgR .* kr) .* besselj(n, wgR .* kc)) ...
        ./ (kr.^2 - kc.^2);
end

