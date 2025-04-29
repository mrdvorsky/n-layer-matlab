clc;
clear;
close all;

%% Inputs
f = @(xs, ys) 1e1 * (1 + xs + 0.1*ys)*0.00 + 1;

nx = 128;
ny = 128;

kRhoLine(:, 1) = linspace(0, 100, 1001);
kPhiLineDeg = 90;

kx(:, 1) = 100 * linspace(-1, 1, 401);
ky(1, :) = 100 * linspace(-1, 1, 401);

Tx(:, 1) = 1 * [0, 1, 0] + 0;
Ty(:, 1) = 1 * [0, 0, 1] + 0;

%% Guassian Grid Points
[nodesX(:, 1), weightsX(:, 1)] = fejer2(nx, 0, 1);
[nodesY(1, :), weightsY(1, :)] = fejer2(ny, 0, 1);

gx = nodesX + 0*nodesY;
gy = nodesY .* (1 - nodesX);
gw = weightsY .* weightsX .* (1 - nodesX);

gx = gx(:);
gy = gy(:);
gw = gw(:);

%% Perform Integrals
samp = f(gx, gy);

spec1 = reshape(...
    nufftn(samp .* gw, [gx, gy]./(2*pi), {kx, ky}), ...
    numel(kx), []);

spec2 = (exp(-1j .* ky) + exp(-1j .* kx) + 1) ./ hypot(kx, ky).^0;

%% Convert Back
[x, y] = fftCoordinates(kx, ky, ApplyFftShift=true);
spat1 = fftshift(ifft2(ifftshift(spec1)));
spat2 = fftshift(ifft2(ifftshift(spec2)));

%% Line Spectrum
specLine = nufftn(samp .* gw, [gx, gy]./(2*pi), ...
    kRhoLine .* [cosd(kPhiLineDeg), sind(kPhiLineDeg)]);

%% Plotting
figure;
showImage(kx, ky, spec1, DisplayFormat="Magnitude");
xlabel("k_x");
ylabel("k_y");

figure;
showImage(kx, ky, spec1 .* hypot(kx, ky), DisplayFormat="MagPhase");
xlabel("k_x");
ylabel("k_y");

figure;
showImage(kx, ky, spec2, DisplayFormat="Magnitude");
xlabel("k_x");
ylabel("k_y");

figure;
showImage(x, y, spat1, DisplayFormat="Magnitude");
xlabel("x");
ylabel("y");

% figure;
% showImage(x, y, spat2, DisplayFormat="Magnitude");
% xlabel("x");
% ylabel("y");

%% Line Plots
figure;
plots(kRhoLine, real(specLine .* kRhoLine.^1), "", LineWidth=1.5);
hold on;
plots(kRhoLine, imag(specLine .* kRhoLine.^1), "", LineWidth=1.5);
grid on;
legend("real", "imag");










