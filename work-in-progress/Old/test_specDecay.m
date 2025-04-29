clc;
clear;
close all;

%% Inputs
x(:, 1) = linspace(-40, 40, 2001);
y(1, :) = linspace(-40, 40, 2001);

r(:, 1) = linspace(0.001, 40, 4001);

nSeries = 15;

%% nLayer
NL = nLayerRectangular(1, 0, waveguideBand="x", convergenceAbsTol=1e-7);

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

%% Calculate Spatial
[kx, ky] = fftCoordinates(x, y);

specEy = specFunEy(kx, ky) + 0*kx.*ky;

scaleFactor = numel(specEy) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
        ./ (2*pi);
Ey = scaleFactor .* ifftshift(ifft2(specEy));

%% Calculate Function
specH = specEy.^2 .* kx.^2;

gamH = scaleFactor .* ifftshift(ifft2(specH));

%% Compute Radial
gamH_interp = griddedInterpolant({x, y}, gamH, "spline");

[phi(1, :), phi_weights(1, :)] = fejer2(1000, 0, 2*pi);
gamH_r = sum(gamH_interp(r .* cos(phi), r .* sin(phi)) .* phi_weights, 2);

%% Compute Multipole Series
R = 2 * hypot(0.5*NL.waveguideA, 0.5*NL.waveguideB);

W = (1 - (r ./ R).^2).^(0:nSeries) .* (r <= R);

c = W \ (gamH_r .* (r <= R));

gamH_fit = W * c;

%% Plot
figure;
showImage(x, y, Ey, DisplayFormat="Real", Normalize=false);
clim([-1, 1] * max(abs(clim())));
colormap colormapPlusMinus;

figure;
showImage(x, y, gamH, DisplayFormat="Imag", Normalize=false);
clim([-1, 1] * max(abs(clim())));
colormap colormapPlusMinus;

figure;
plot(r, db(gamH_r), "", LineWidth=1.5);
grid on;

figure;
plot(r, abs(gamH_r), "", LineWidth=1.5);
hold on;
plot(r, abs(gamH_fit), ":", LineWidth=1.5);
grid on;




