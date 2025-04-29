clc;
clear;
close all;

%% Inputs
wgR = 1.0;

x(:, 1) = 20.5 * linspace(-1, 1, 1001);
y(1, :) = 20.5 * linspace(-1, 1, 1001);

phiSample(:, 1) = linspace(0, 2*pi, 1001);

m = 1;
n = 0;

kc = besselj_zeros(n, m) ./ wgR;
% kc = besselj_der_zeros(n, m) ./ wgR;

%% Calculate Spectrums
[kx, ky] = fftCoordinates(x, y, ApplyFftShift=true);
[xs, ys, weights] = gaussCircle(wgR, 401, 401);

ExSamp = fieldEx(xs, ys, wgR, kc, n);
EySamp = fieldEy(xs, ys, wgR, kc, n);
EzSamp = fieldEz(xs, ys, wgR, kc, n);

ExSpec = reshape(nufftn(weights .* ExSamp, [xs, ys] ./ (2*pi), {kx, ky}), numel(kx), []);
EySpec = reshape(nufftn(weights .* EySamp, [xs, ys] ./ (2*pi), {kx, ky}), numel(kx), []);
EzSpec = reshape(nufftn(weights .* EzSamp, [xs, ys] ./ (2*pi), {kx, ky}), numel(kx), []);

%% Analytical Spectrums
kphi = angle(kx + 1j*ky);
kr = hypot(kx, ky);

EzSpecAn = 2 * cos(n .* kphi) ...
        .* (kc .* besselj(n, wgR .* kr) .* besselj(n - 1, wgR .* kc) ...
        - kr .* besselj(n - 1, wgR .* kr) .* besselj(n, wgR .* kc)) ...
        ./ (kr.^2 - kc.^2);

ExSpecAn = -ky .* EzSpecAn;
EySpecAn = kx .* EzSpecAn;

ExSpecDiff2 = (sin((n - 1) .* kphi) .* besselj(n - 1, wgR .* kr) ...
    + sin((n + 1) .* kphi) .* besselj(n + 1, wgR .* kr)) ...
    .* (besselj(n, kc .* wgR));
EySpecDiff2 = (cos((n - 1) .* kphi) .* besselj(n - 1, wgR .* kr) ...
    - cos((n + 1) .* kphi) .* besselj(n + 1, wgR .* kr)) ...
    .* (besselj(n, kc .* wgR));


% ExSpecAn = 1j .* ky .* EzSpec;
% EySpecAn = -1j .* kx .* EzSpec;

ExSpec2 = ExSpecAn + ExSpecDiff2;
EySpec2 = EySpecAn + EySpecDiff2;


sm = sin((n - 1).*kphi);
sp = sin((n + 1).*kphi);
cm = cos((n - 1).*kphi);
cp = cos((n + 1).*kphi);
Jnkr = besselj(n, wgR .* kr);
Jnkc = besselj(n, wgR .* kc);
Jmkr = besselj(n - 1, wgR .* kr);
Jmkc = besselj(n - 1, wgR .* kc);
Jpkr = besselj(n + 1, wgR .* kr);
Jpkc = besselj(n + 1, wgR .* kc);
JnkrOverKr = Jnkr ./ kr;
JnkrOverKr(kr == 0) = 0.5*wgR*(n == 1);

% ExSpec2 = (kr.*kc.*(sm - sp) .* Jnkr.*Jmkc ...
%     - kc.^2 .* Jnkc .* (sm.*Jmkr + sp.*Jpkr) ...
%     + sp.*(2*n/wgR).*kr .* Jnkr.*Jnkc) ...
%     ./ (kr.^2 - kc.^2);

ExSpec2 = kc .* (sm - sp) .* (kr.*Jnkr.*Jmkc - kc.*Jmkr.*Jnkc) ...
    ./ (kr.^2 - kc.^2) ...
    + (2*n./wgR) .* sp .* JnkrOverKr.*Jnkc;

EySpec2 = kc .* (cm + cp) .* (kr.*Jnkr.*Jmkc - kc.*Jmkr.*Jnkc) ...
    ./ (kr.^2 - kc.^2) ...
    - (2*n./wgR) .* cp .* JnkrOverKr.*Jnkc;

scaleFac = (-1j).^(n + 1) * pi .* wgR;
ExSpec2 = ExSpec2 .* scaleFac;
EySpec2 = EySpec2 .* scaleFac;

scaleX = ExSpec(:) \ ExSpec2(:)
scaleY = EySpec(:) \ EySpec2(:)

%% Magnitudes
scaleKxKy = abs((kx(1) - kx(2)) .* (ky(1) - ky(2)));
modeMag = sum(abs(ExSpec2).^2 + abs(EySpec2).^2, "all") .* scaleKxKy

magTest = 8*pi.^2 * pi * integral(@(x) x .* kc.^2 .* besselj_der(n, kc.*x).^2 ...
    + n.^2 .* besselj(n, kc.*x).^2 ./ x, ...
    0, wgR)

%% Plot Spectrums
dispFormat = "Magnitude";
dispFormat = "MagPhase";

figure;
showImage(kx, ky, ExSpec, DisplayFormat=dispFormat);
title("Ex");

figure;
showImage(kx, ky, ExSpec2, DisplayFormat=dispFormat);
title("Ex'");


figure;
showImage(kx, ky, EySpec, DisplayFormat=dispFormat);
title("Ey");

figure;
showImage(kx, ky, EySpec2, DisplayFormat=dispFormat);
title("Ey'");









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
    y = (-1j).^(n) * 2 * pi * wgR ...
        .* (kc .* besselj(n, wgR .* kr) .* besselj(n - 1, wgR .* kc) ...
        - kr .* besselj(n - 1, wgR .* kr) .* besselj(n, wgR .* kc)) ...
        ./ (kr.^2 - kc.^2);
end

function [Ex] = fieldEx(x, y, wgR, kc, n)
    r = hypot(x, y);
    phi = angle(x + 1j*y);

    modeWin = (r <= wgR);
    Ephi = -kc .* besselj_der(n, kc * r) .* cos(n .* phi) .* modeWin;
    Er = -n .* besselj(n, kc * r) .* sin(n .* phi) ./ r .* modeWin;
    Er(r == 0) = 0.5*kc*(n == 1);
    
    Ex = Er .* cos(phi) - Ephi .* sin(phi);
end

function [Ey] = fieldEy(x, y, wgR, kc, n)
    r = hypot(x, y);
    phi = angle(x + 1j*y);

    modeWin = (r <= wgR);
    Ephi = -kc .* besselj_der(n, kc * r) .* cos(n .* phi) .* modeWin;
    Er = -n .* besselj(n, kc * r) .* sin(n .* phi) ./ r .* modeWin;
    Er(r == 0) = 0.5*kc*(n == 1);
    
    Ey = Er .* sin(phi) + Ephi .* cos(phi);
end

function [Ez] = fieldEz(x, y, wgR, kc, n)
    r = hypot(x, y);
    phi = angle(x + 1j*y);

    modeWin = (r <= wgR);
    Ez = besselj(n, kc * r) .* cos(n .* phi) .* modeWin;
end
