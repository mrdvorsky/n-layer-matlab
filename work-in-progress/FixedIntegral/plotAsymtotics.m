clc;
clear;
close all;

%% Inputs
krMax = 1000;
Nkr = 500;

Nphi = 10;

%% nLayer
wgA = 1;
wgB = 1/2;

%% Get Spectrum Functions
[spEx1, spEy1] = nLayerOpenEnded.getSpectrumRectangular(...
    wgA, wgB, 3, 2, "TM");
[spEx2, spEy2] = nLayerOpenEnded.getSpectrumRectangular(...
    wgA, wgB, 3, 0, "TM");

kr(:, 1) = linspace(0.5 * krMax, krMax, Nkr);
phi(1, :) = fejer2(Nphi, 0, pi/2);
kx = kr .* cos(phi);
ky = kr .* sin(phi);

spec = spEy1(kx, ky);

%% Fit
pointsX(1, 1, :) = 0.5 * [-1, 1, -1, 1, -1, 0, 0, 1];
pointsY(1, 1, :) = 0.25 * [1, 1, -1, -1, 0, 1, -1, 0];
% pointsX(1, 1, :) = 0.5 * [-1, 1, -1, 1];
% pointsY(1, 1, :) = 0.25 * [1, 1, -1, -1];
% krPow(1, 1, 1, :) = [1, 2];
kxPow(1, 1, 1, :) = [1, 2, 2];
kyPow(1, 1, 1, :) = [2, 1, 2];
% kxPow(1, 1, 1, :) = [1, 1, 1, 2, 2, 2, 3, 3, 3];
% kyPow(1, 1, 1, :) = [1, 2, 3, 1, 2, 3, 1, 2, 3];

% pointsX(1, 1, :) = 0.5 * [-1, 0, 0, 1];
% pointsY(1, 1, :) = 0.25 * [0, 1, -1, 0];
% kxPow(1, 1, 1, :) = [1, 1, 2, 2];
% kyPow(1, 1, 1, :) = [1, 2, 1, 2];

f = exp(-1j .* kx .* pointsX) .* exp(-1j .* ky .* pointsY) ./ kx.^(kxPow) ./ ky.^(kyPow);
f = f ./ rms(f, [1, 2]) .* max(abs(spec(:)));
% f = exp(-1j .* kx .* pointsX) .* exp(-1j .* ky .* pointsY) ./ hypot(kx, ky).^krPow;
% f = exp(-1j .* hypot(kx, ky) .* hypot(pointsY, pointsX)) ./ hypot(kx, ky).^krPow;

A = reshape(f, numel(spec), []);

coeff = A \ spec(:)
db(coeff)

spec_fit = reshape(A * coeff, size(spec));

%% Sample
figure;
semilogx(kr, db(spec(:, 10)), "", LineWidth=1.5);
hold on;
semilogx(kr, db(spec_fit(:, 10)), "x", LineWidth=1.5);
grid on;

figure;
semilogx(kr, db(spec(:, :) - spec_fit(:, :)) - db(max(abs(spec), [], "all")), "", LineWidth=1.5);
grid on;









