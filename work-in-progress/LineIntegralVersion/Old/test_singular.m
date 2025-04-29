clc;
clear;
close all;

%% Inputs
% x(:, 1) = linspace(-0.5, 0.5, 11);

L = 1;
Nkr = 16000;
Nphi = 2000;

xL = 1;

%% Integration 1
[kr(1, :, 1), kr_w(1, :, 1)] = fejer2_halfOpen(Nkr, L);
[kphi(1, 1, :), kphi_w(1, 1, :)] = fejer2(Nphi, 0, pi/2);

I1 = sum(sinc(xL * kr .* cos(kphi) ./ pi) ./ (1 + kr) ...
    .* kphi_w .* kr_w, [2, 3]) ./ pi

%% Integration 2
[x(:, 1), x_w(:, 1)] = fejer2(100, -1, 1);

sum();


%% Plotting
xp(:, 1) = linspace(-10, 10, 4001);
yp(1, :) = linspace(-10, 10, 4001);
[kx_p, ky_p] = fftCoordinates(xp, yp);


spec = sinc(kx_p./pi).^1 ./ (1 + hypot(kx_p, ky_p));
% spec = sinc(kx_p./pi);


scale = abs(kx_p(2) - kx_p(1)) * abs(ky_p(2) - ky_p(1));
spat1 = scale * fftshift(fft2(spec + 0*kx_p*ky_p));
figure;
showImage(xp, yp, spat1, DisplayFormat="Magnitude");



