clc;
clear;
close all;

%% Inputs
wgR = 1;

x(:, 1) = 2.0 * linspace(-1, 1, 4001);
y(1, :) = 2.0 * linspace(-1, 1, 4001);

m = 2;
n = 3;

kc = besselj_der_zeros(n, m) ./ wgR;

%% Coordinates
r = hypot(x, y);
phi = angle(x + 1j*y);

dx = abs(x(1) - x(2));
dy = abs(y(1) - y(2));

%% Calculate Fields
modeWin = (r <= wgR);

Ephi = kc .* besselj_der(n, kc * r) .* cos(n .* phi) .* modeWin;
Er = n .* besselj(n, kc * r) .* sin(n .* phi) ./ r .* modeWin;
Er(r == 0) = 0.5*(n == 1);

Ex = Er .* cos(phi) - Ephi .* sin(phi);
Ey = Er .* sin(phi) + Ephi .* cos(phi);

%% Spectrums
[kx, ky] = fftCoordinates(x, y, ApplyFftShift=true);
kr = hypot(kx, ky);
kphi = angle(kx + 1j*ky);

scaleFactor = abs(x(1) - x(2)) * abs(y(1) - y(2));
Ex_hat = -1j*fftshift(fft2(ifftshift(Ex))) .* scaleFactor;
Ey_hat = -1j*fftshift(fft2(ifftshift(Ey))) .* scaleFactor;

% Ex_der_hat = -1j*fftshift(fft2(ifftshift(Ex_der))) .* scaleFactor;
% Ey_der_hat = -1j*fftshift(fft2(ifftshift(Ey_der))) .* scaleFactor;
[phiInt, weights_phi] = fejer2(1000, 0, 0.5*pi);
phiInt = [phiInt; phiInt + 0.5*pi; phiInt + 1.0*pi; phiInt + 1.5*pi];
weights_phi = repmat(weights_phi, 4, 1);


xInt = wgR .* cos(phiInt);
yInt = wgR .* sin(phiInt);
Ex_derInt = besselj(n, kc * wgR) .* cos(n * phiInt) .* sign(-sin(phiInt));
Ey_derInt = besselj(n, kc * wgR) .* cos(n * phiInt) .* sign(-cos(phiInt));
Ex_der_hat = reshape(nufftn(weights_phi .* Ex_derInt, [xInt, yInt]./(2*pi), {kx, ky}), ...
    numel(kx), numel(ky)) .* (-1j*scaleFactor);
Ey_der_hat = reshape(nufftn(weights_phi .* Ey_derInt, [xInt, yInt]./(2*pi), {kx, ky}), ...
    numel(kx), numel(ky)) .* (-1j*scaleFactor);

Ex_hat2 = (2*pi) * kr .* sin(kphi) .* cos(n * kphi) .* besselInt(kr, n, kc, wgR);
Ey_hat2 = -(2*pi) * kr .* cos(kphi) .* cos(n * kphi) .* besselInt(kr, n, kc, wgR);

errX_rel = abs(max(Ex_hat - Ex_hat2, [], "all") ./ max(Ex_hat2, [], "all"))
errY_rel = abs(max(Ey_hat - Ey_hat2, [], "all") ./ max(Ey_hat2, [], "all"))

%% Diff
Ex_der_hat2 = Ex_hat - Ex_hat2;
Ey_der_hat2 = Ey_hat - Ey_hat2;

Ex_der = fftshift(ifft2(ifftshift(Ex_der_hat2)));
Ey_der = fftshift(ifft2(ifftshift(Ey_der_hat2)));

%% Plot 1D
Ex_der_grid = griddedInterpolant({x, y}, Ex_der, "spline");
Ey_der_grid = griddedInterpolant({x, y}, Ey_der, "spline");

phiSample = linspace(0, 2*pi, 1001);
xS = wgR * cos(phiSample);
yS = wgR * sin(phiSample);

ExLine = Ex_der_grid(xS, yS);
EyLine = Ey_der_grid(xS, yS);

figure;
plot(rad2deg(phiSample), real(ExLine) + imag(ExLine), "", Linewidth=1.5);
hold on;
plot(rad2deg(phiSample), real(EyLine) + imag(EyLine), "", Linewidth=1.5);


%% Plotting
scale = max(abs(cat(3, Ex, Ey)), [], "all");
dispForm = "Magnitude";

figure;
subplot(2, 2, 1);
showImage(x, y, Ex, DisplayFormat=dispForm);
clim(scale * [0, 1]);
title("E_{x}");

subplot(2, 2, 2);
showImage(x, y, Ey, DisplayFormat=dispForm);
clim(scale * [0, 1]);
title("E_{y}");

scale2 = max(abs(cat(3, Ex_der, Ey_der)), [], "all");
subplot(2, 2, 3);
showImage(x, y, abs(Ex_der), DisplayFormat=dispForm);
clim(scale2 * [0, 1]);
title("E_{x}");

subplot(2, 2, 4);
showImage(x, y, abs(Ex_der), DisplayFormat=dispForm);
clim(scale2 * [0, 1]);
title("E_{y}");

%% Plotting
scale = max(abs(cat(3, Ex_hat, Ey_hat, Ex_hat2, Ey_hat2)), [], "all");
dispForm = "Magnitude";

figure;
subplot(2, 2, 1);
showImage(kx, ky, Ex_hat, DisplayFormat=dispForm);
clim(scale * [0, 1]);
title("E_{x}");

subplot(2, 2, 2);
showImage(kx, ky, Ey_hat, DisplayFormat=dispForm);
clim(scale * [0, 1]);
title("E_{y}");

subplot(2, 2, 3);
showImage(kx, ky, Ex_hat2, DisplayFormat=dispForm);
clim(scale * [0, 1]);
title("E_{x}");

subplot(2, 2, 4);
showImage(kx, ky, Ey_hat2, DisplayFormat=dispForm);
clim(scale * [0, 1]);
title("E_{y}");

%% Plot Difference
Ex_hat_diff = Ex_hat - Ex_hat2;
Ey_hat_diff = Ey_hat - Ey_hat2;

scale = max(abs(cat(3, Ex_hat_diff, Ey_hat_diff)), [], "all");
dispForm = "Magnitude";

figure;
subplot(2, 2, 1);
showImage(kx, ky, Ex_hat_diff ./ 1, DisplayFormat=dispForm);
clim(scale * [0, 1]);
title("E_{x}");

subplot(2, 2, 2);
showImage(kx, ky, Ey_hat_diff ./ 1, DisplayFormat=dispForm);
clim(scale * [0, 1]);
title("E_{y}");

scale2 = max(abs(cat(3, Ex_der_hat, Ey_der_hat)), [], "all");
subplot(2, 2, 3);
showImage(kx, ky, Ex_der_hat ./ 1, DisplayFormat=dispForm);
clim(scale2 * [0, 1]);
title("E_{x}");

subplot(2, 2, 4);
showImage(kx, ky, Ey_der_hat ./ 1, DisplayFormat=dispForm);
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

function [Ex] = fieldEx(x, y, wgR, kc, n)
    r = hypot(x, y);
    phi = angle(x + 1j*y);

    modeWin = (r <= wgR);
    Ephi = kc .* besselj_der(n, kc * r) .* cos(n .* phi) .* modeWin;
    Er = n .* besselj(n, kc * r) .* sin(n .* phi) ./ r .* modeWin;
    Er(r == 0) = 0.5*(n == 1);
    
    Ex = Er .* cos(phi) - Ephi .* sin(phi);
end

function [Ey] = fieldEy(x, y, wgR, kc, n)
    r = hypot(x, y);
    phi = angle(x + 1j*y);

    modeWin = (r <= wgR);
    Ephi = kc .* besselj_der(n, kc * r) .* cos(n .* phi) .* modeWin;
    Er = n .* besselj(n, kc * r) .* sin(n .* phi) ./ r .* modeWin;
    Er(r == 0) = 0.5*(n == 1);
    
    Ey = Er .* sin(phi) + Ephi .* cos(phi);
end

