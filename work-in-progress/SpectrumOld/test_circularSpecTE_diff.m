clc;
clear;
close all;

%% Inputs
wgR = 1;

x(:, 1) = 20.5 * linspace(-1, 1, 1001);
y(1, :) = 20.5 * linspace(-1, 1, 1001);

phiSample(:, 1) = linspace(0, 2*pi, 1001);

m = 1;
n = 0;

kc = besselj_der_zeros(n, m) ./ wgR;

%% Calculate Spectrums
[kx, ky] = fftCoordinates(x, y, ApplyFftShift=true);

[xs, ys, weights] = gaussCircle(wgR, 401, 401);

ExSamp = fieldEx(xs, ys, wgR, kc, n);
EySamp = fieldEy(xs, ys, wgR, kc, n);
EzSamp = fieldEz(xs, ys, wgR, kc, n);

ExSpec = reshape(nufftn(weights .* ExSamp, [xs, ys] ./ (2*pi), {kx, ky}), numel(kx), []);
EySpec = reshape(nufftn(weights .* EySamp, [xs, ys] ./ (2*pi), {kx, ky}), numel(kx), []);
EzSpec = reshape(nufftn(weights .* EzSamp, [xs, ys] ./ (2*pi), {kx, ky}), numel(kx), []);

ExSpec2 = -1j*ky .* EzSpec;
EySpec2 = 1j*kx .* EzSpec;

%% Differences
xSample = wgR .* cos(phiSample);
ySample = wgR .* sin(phiSample);

ExSpecDiff = ExSpec - ExSpec2;
EySpecDiff = EySpec - EySpec2;

kr = hypot(kx, ky);
weightFun = kr < sqrt(0.5)*max(kr, [], "all");
scaleKxKy = 1 ./ abs((kx(1) - kx(2)) .* (ky(1) - ky(2)));
ExDiff = nufftn(ExSpecDiff .* weightFun, {kx, ky}, -[xSample, ySample]./(2*pi)) .* scaleKxKy;
EyDiff = nufftn(EySpecDiff .* weightFun, {kx, ky}, -[xSample, ySample]./(2*pi)) .* scaleKxKy;

fit_ExDiff = real(ExDiff);
fit_EyDiff = real(EyDiff);
fit_phi = phiSample;

%% Analytical Diff
kphi = angle(kx + 1j*ky);
kr = hypot(kx, ky);

ExSpecDiff2 = (-sin((n - 1) .* kphi) .* besselj(n - 1, wgR .* kr) ...
    - sin((n + 1) .* kphi) .* besselj(n + 1, wgR .* kr)) ...
    .* (besselj(n, kc .* wgR) .* wgR .* pi);
EySpecDiff2 = (cos((n - 1) .* kphi) .* besselj(n - 1, wgR .* kr) ...
    - cos((n + 1) .* kphi) .* besselj(n + 1, wgR .* kr)) ...
    .* (besselj(n, kc .* wgR) .* wgR .* pi);

scaleX = max(ExSpecDiff, [], "all") ./ max(ExSpecDiff2, [], "all")
scaleY = max(EySpecDiff, [], "all") ./ max(EySpecDiff2, [], "all")

%% Plot Spectrums
dispFormat = "Magnitude";
% dispFormat = "MagPhase";

figure;
showImage(kx, ky, ExSpec, DisplayFormat=dispFormat);
title("Ex");

figure;
showImage(kx, ky, EySpec, DisplayFormat=dispFormat);
title("Ey");

figure;
showImage(kx, ky, EzSpec, DisplayFormat=dispFormat);
title("Ez");

figure;
showImage(kx, ky, ExSpecDiff, DisplayFormat=dispFormat);
title("Ex'");

figure;
showImage(kx, ky, EySpecDiff, DisplayFormat=dispFormat);
title("Ey'");

figure;
showImage(kx, ky, ExSpecDiff2, DisplayFormat=dispFormat);
title("Ex'");

figure;
showImage(kx, ky, EySpecDiff2, DisplayFormat=dispFormat);
title("Ey'");

%% Plot Spatial
figure;
plot(rad2deg(phiSample), real(ExDiff), "", LineWidth=1.5);
hold on;
plot(rad2deg(phiSample), imag(ExDiff), "", LineWidth=1.5);
grid on;
legend("Real", "Imag");
title("ExDiff");

figure;
plot(rad2deg(phiSample), real(EyDiff), "", LineWidth=1.5);
hold on;
plot(rad2deg(phiSample), imag(EyDiff), "", LineWidth=1.5);
grid on;
legend("Real", "Imag");
title("EyDiff");

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

function [Ez] = fieldEz(x, y, wgR, kc, n)
    r = hypot(x, y);
    phi = angle(x + 1j*y);

    modeWin = (r <= wgR);
    Ez = besselj(n, kc * r) .* cos(n .* phi) .* modeWin;
end
