clc;
clear;
close all;

%% Inputs
wgR = 1;

x(:, 1) = 50.1 * linspace(-1, 1, 2001);
y(1, :) = 50.1 * linspace(-1, 1, 2001);

m = 1;
n = 1;

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

Ex_hat = 1j*fftshift(fft2(ifftshift(Ex)));
Ey_hat = 1j*fftshift(fft2(ifftshift(Ey)));

Er_hat = 1j*fftshift(fft2(ifftshift(Er)));
Ephi_hat = 1j*fftshift(fft2(ifftshift(Ephi)));

Ex_hat2 = Er_hat .* cos(kphi) - Ephi_hat .* sin(kphi);
Ey_hat2 = Er_hat .* sin(kphi) + Ephi_hat .* cos(kphi);

%% Plotting
scale = max(abs(cat(3, Ex, Ey, Er, Ephi)), [], "all");

figure;
subplot(2, 2, 1);
showImage(x, y, Er);
clim(scale * [-1, 1]);
title("E_{r}");

subplot(2, 2, 2);
showImage(x, y, Ephi);
clim(scale * [-1, 1]);
title("E_{\phi}");

subplot(2, 2, 3);
showImage(x, y, Ex);
clim(scale * [-1, 1]);
title("E_{x}");

subplot(2, 2, 4);
showImage(x, y, Ey);
clim(scale * [-1, 1]);
title("E_{y}");

%% Plotting
scale = max(abs(cat(3, Ex_hat, Ey_hat, Er_hat, Ephi_hat)), [], "all");

figure;
subplot(3, 2, 1);
showImage(kx, ky, Er_hat, DisplayFormat="Magnitude");
clim(scale * [0, 1]);
title("E_{r}");

subplot(3, 2, 2);
showImage(kx, ky, Ephi_hat, DisplayFormat="Magnitude");
clim(scale * [0, 1]);
title("E_{\phi}");

subplot(3, 2, 3);
showImage(kx, ky, Ex_hat, DisplayFormat="Magnitude");
clim(scale * [0, 1]);
title("E_{x}");

subplot(3, 2, 4);
showImage(kx, ky, Ey_hat, DisplayFormat="Magnitude");
clim(scale * [0, 1]);
title("E_{y}");

subplot(3, 2, 5);
showImage(kx, ky, Ex_hat2, DisplayFormat="Magnitude");
clim(scale * [0, 1]);
title("E_{x}");

subplot(3, 2, 6);
showImage(kx, ky, Ey_hat2, DisplayFormat="Magnitude");
clim(scale * [0, 1]);
title("E_{y}");

%% Plotting
figure;
plot(kx, abs(Ex_hat(:, ky == 0)), "", LineWidth=1.5);
hold on;
plot(kx, abs(Ex_hat2(:, ky == 0)), "", LineWidth=1.5);
grid on;
xlabel("k_x");
ylabel("E_{r}");

%% Arrow Plots
% step = 30;
% figure;
% quiver(x(1:step:end, 1:step:end), y(1:step:end, 1:step:end), ...
%     Ex(1:step:end, 1:step:end).', Ey(1:step:end, 1:step:end).');
% axis image;
% grid on;


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

