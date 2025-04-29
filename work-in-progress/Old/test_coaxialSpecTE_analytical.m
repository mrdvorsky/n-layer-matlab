clc;
clear;
close all;

%% Inputs
wgR_inner = 0.5;
wgR_outer = 1.0;

x(:, 1) = 5.5 * linspace(-1, 1, 1001);
y(1, :) = 5.5 * linspace(-1, 1, 1001);

phiSample(:, 1) = linspace(0, 2*pi, 1001);

m = 1;
n = 1;

[kc, a, b] = coaxial_cutoffs_TM(wgR_inner, wgR_outer, m, n);
% [kc, a, b] = coaxial_cutoffs_TE(wgR_inner, wgR_outer, m, n);

%% Calculate Spectrums
[kx, ky] = fftCoordinates(x, y, ApplyFftShift=true);
[xs, ys, weights] = gaussCircle(wgR_outer, 401, 401, wgR_inner);

ExSamp = fieldEx(xs, ys, wgR_inner, wgR_outer, kc, a, b, m);
EySamp = fieldEy(xs, ys, wgR_inner, wgR_outer, kc, a, b, m);
EzSamp = fieldEz(xs, ys, wgR_inner, wgR_outer, kc, a, b, m);

ExSpec = reshape(nufftn(weights .* ExSamp, [xs, ys] ./ (2*pi), {kx, ky}), numel(kx), []);
EySpec = reshape(nufftn(weights .* EySamp, [xs, ys] ./ (2*pi), {kx, ky}), numel(kx), []);
EzSpec = reshape(nufftn(weights .* EzSamp, [xs, ys] ./ (2*pi), {kx, ky}), numel(kx), []);

%% Analytical Spectrums
kphi = angle(kx + 1j*ky);
kr = hypot(kx, ky);

EzSpecAn = 2 * cos(m .* kphi) ...
        .* (kc .* besselj(m, wgR_outer .* kr) .* besselj(m - 1, wgR_outer .* kc) ...
        - kr .* besselj(m - 1, wgR_outer .* kr) .* besselj(m, wgR_outer .* kc)) ...
        ./ (kr.^2 - kc.^2);

% ExSpecAn = -ky .* EzSpecAn;
% EySpecAn = kx .* EzSpecAn;

ExSpecAn = kx .* EzSpecAn;
EySpecAn = ky .* EzSpecAn;

ExSpecDiff2 = (sin((m - 1) .* kphi) .* besselj(m - 1, wgR_outer .* kr) ...
    + sin((m + 1) .* kphi) .* besselj(m + 1, wgR_outer .* kr)) ...
    .* (besselj(m, kc .* wgR_outer));
EySpecDiff2 = (cos((m - 1) .* kphi) .* besselj(m - 1, wgR_outer .* kr) ...
    - cos((m + 1) .* kphi) .* besselj(m + 1, wgR_outer .* kr)) ...
    .* (besselj(m, kc .* wgR_outer));


ExSpecAn = 1j .* kx .* EzSpec;
EySpecAn = -1j .* ky .* EzSpec;

ExSpec2 = ExSpecAn;
EySpec2 = EySpecAn;


sm = sin((m - 1).*kphi);
sp = sin((m + 1).*kphi);
cm = cos((m - 1).*kphi);
cp = cos((m + 1).*kphi);
Jnkr = besselj(m, wgR_outer .* kr);
Jnkc = besselj(m, wgR_outer .* kc);
Jmkr = besselj(m - 1, wgR_outer .* kr);
Jmkc = besselj(m - 1, wgR_outer .* kc);
Jpkr = besselj(m + 1, wgR_outer .* kr);
Jpkc = besselj(m + 1, wgR_outer .* kc);
JnkrOverKr = Jnkr ./ kr;
JnkrOverKr(kr == 0) = 0.5*wgR_outer*(m == 1);

% ExSpec2 = kc .* (sm - sp) .* (kr.*Jnkr.*Jmkc - kc.*Jmkr.*Jnkc) ...
%     ./ (kr.^2 - kc.^2) ...
%     + (2*m./wgR_outer) .* sp .* JnkrOverKr.*Jnkc;
% 
% EySpec2 = kc .* (cm + cp) .* (kr.*Jnkr.*Jmkc - kc.*Jmkr.*Jnkc) ...
%     ./ (kr.^2 - kc.^2) ...
%     - (2*m./wgR_outer) .* cp .* JnkrOverKr.*Jnkc;

scaleFac = (-1j).^(m + 1) * pi .* wgR_outer;
ExSpec2 = ExSpec2 .* scaleFac;
EySpec2 = EySpec2 .* scaleFac;

scaleX = ExSpec(:) \ ExSpec2(:)
scaleY = EySpec(:) \ EySpec2(:)

%% Plot Spectrums
dispFormat = "Magnitude";
% dispFormat = "MagPhase";

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

%% Plot Fields
dispFormat = "Magnitude";
% dispFormat = "MagPhase";

[x, y] = fftCoordinates(kx, ky, ApplyFftShift=true);
phi = angle(x + 1j*y);

Ex = fieldEx(x, y, wgR_inner, wgR_outer, kc, a, b, m);
Ey = fieldEy(x, y, wgR_inner, wgR_outer, kc, a, b, m);
Ez = fieldEz(x, y, wgR_inner, wgR_outer, kc, a, b, m);

figure;
showImage(x, y, Ez, DisplayFormat=dispFormat);
title("Ez");











%% Field Definitions
function [Ex] = fieldEx(x, y, wgR1, wgR2, kc, a, b, m)
    r = hypot(x, y);
    phi = angle(x + 1j*y);

    modeWin = (r <= wgR2) & (r >= wgR1);
    Er = -kc .* besseljy_derivative(a, b, m, kc * r) .* cos(m .* phi) .* modeWin;
    Ephi = -m .* besseljy(a, b, m, kc * r) .* sin(m .* phi) ./ r .* modeWin;
    
    Ex = Er .* cos(phi) - Ephi .* sin(phi);
end

function [Ey] = fieldEy(x, y, wgR1, wgR2, kc, a, b, m)
    r = hypot(x, y);
    phi = angle(x + 1j*y);

    modeWin = (r <= wgR2) & (r >= wgR1);
    Er = -kc .* besseljy_derivative(a, b, m, kc * r) .* cos(m .* phi) .* modeWin;
    Ephi = -m .* besseljy(a, b, m, kc * r) .* sin(m .* phi) ./ r .* modeWin;
    
    Ey = Er .* sin(phi) + Ephi .* cos(phi);
end

function [Ez] = fieldEz(x, y, wgR1, wgR2, kc, a, b, m)
    r = hypot(x, y);
    phi = angle(x + 1j*y);

    modeWin = (r <= wgR2) & (r >= wgR1);
    Ez = besseljy(a, b, m, kc * r) .* cos(m .* phi) .* modeWin;
end

%% Bessel Functions
function [y] = besseljy(a, b, v, x)
    y = a.*besselj(v, x) + b.*bessely(v, x);
end

function [y] = besseljy_derivative(a, b, v, x)
    y = 0.5 * (besseljy(a, b, v - 1, x) - besseljy(a, b, v + 1, x));
end

function [kc, a, b] = coaxial_cutoffs_TM(wgR1, wgR2, m, n)
    kc = besseljy_det_zeros(m, n, wgR1 ./ wgR2) ./ wgR2;
    ab = null([besseljy(1, 0, m, kc .* wgR2), ...
        besseljy(0, 1, m, kc .* wgR2)]);
    a = ab(1);
    b = ab(2);
end

function [kc, a, b] = coaxial_cutoffs_TE(wgR1, wgR2, m, n)
    kc = besseljy_det_derivative_zeros(m, n, wgR1 ./ wgR2) ./ wgR2;
    ab = null([besseljy_derivative(1, 0, m, kc .* wgR2), ...
        besseljy_derivative(0, 1, m, kc .* wgR2)]);
    a = ab(1);
    b = ab(2);
end


%% Bessel Function Zeros
function [x] = besseljy_det_zeros(v, n, a)
    fun = @(y) besseljy_det(v, abs(y), a);
    x = abs(fzero(fun, (v + (v == 0))));
    for ii = 2:n
        x_guess = x + (0.5*pi ./ (1 - a))*[0.9, 1.1];
        while sum(sign(fun(x_guess))) ~= 0
            x_guess(2) = x + (x_guess(2) - x)*1.1;
        end
        x = fzero(fun, x_guess);
    end
end

function [x] = besseljy_det_derivative_zeros(v, n, a)
    fun = @(y) besseljy_derivative_det(v, abs(y), a);
    x = abs(fzero(fun, (v + (v == 0))));
    for ii = 2:n
        x_guess = x + (0.5*pi ./ (1 - a))*[0.9, 1.1];
        while sum(sign(fun(x_guess))) ~= 0
            x_guess(2) = x + (x_guess(2) - x)*1.1;
        end
        x = fzero(fun, x_guess);
    end
end

function [y] = besseljy_det(v, x, a)
    y = besselj(v, x).*bessely(v, a.*x) ...
        - besselj(v, a.*x).*bessely(v, x);
end

function [y] = besseljy_derivative_det(v, x, a)
    y = 0.5 * (besseljy_det(v - 1, x, a) - besseljy_det(v + 1, x, a));
end


