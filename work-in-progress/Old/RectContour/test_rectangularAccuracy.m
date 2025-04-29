clc;
clear;
close all;

%% Inputs
wgA = 1;
wgB = 0.5;

m = 0;
n = 1;

kx(:, 1) = 400 * linspace(-1, 1, 1001);
ky(1, :) = 400 * linspace(-1, 1, 1001);

nx = 800;
ny = 800;

%% Calculate Fourier Transform Analytically
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(wgA, wgB, m, n, "TE");

specEx = specFunEx(kx, ky) + 0*kx.*ky;
specEy = specFunEy(kx, ky) + 0*kx.*ky;

%% Calculate Numerically
ExFun = @(x, y) -1j./(2*pi) .* m ./ wgB .* sqrt(4 ./ (wgA .* wgB)) ./ hypot(m ./ wgB, n ./ wgA) ...
    .* cos(m*pi./wgA .* (x - 0.5*wgA)) ...
    .* sin(n*pi./wgB .* (y - 0.5*wgB));

EyFun = @(x, y) 1j./(2*pi) .* n ./ wgA .* sqrt(4 ./ (wgA .* wgB)) ./ hypot(m ./ wgB, n ./ wgA) ...
    .* sin(m*pi./wgA .* (x - 0.5*wgA)) ...
    .* cos(n*pi./wgB .* (y - 0.5*wgB));

[xs(:, 1), xs_weights(:, 1)] = fejer2(nx, -0.5*wgA, 0.5*wgA);
[ys(1, :), ys_weights(1, :)] = fejer2(ny, -0.5*wgB, 0.5*wgB);

ExSamp = ExFun(xs, ys) .* xs_weights .* ys_weights;
tic;
specEx_num = reshape(...
    nufftn(ExSamp, {xs ./ (2*pi), ys ./ (2*pi)}, {kx, ky}), ...
    numel(kx), numel(ky));
toc;

EySamp = EyFun(xs, ys) .* xs_weights .* ys_weights;
tic;
specEy_num = reshape(...
    nufftn(EySamp, {xs ./ (2*pi), ys ./ (2*pi)}, {kx, ky}), ...
    numel(kx), numel(ky));
toc;

%% Calculate Conv
kPhi = angle(kx + 1j*ky);
convSpec = (specEx .* cos(kPhi) + specEy .* sin(kPhi)).^2;
convSpec = (specEx .* ky + specEy .* kx).^2;

gFilt = 1 + 0*exp(-hypot(kx, ky).^2 ./ (100.^2));

[x, y] = fftCoordinates(kx, ky, ApplyFftShift=true);
convVal = ifftshift(ifft2(ifftshift(specEy.^2 .* kx.^2 .* gFilt)));
convVal = ifftshift(ifft2(ifftshift(convSpec .* gFilt)));

%% Plotting
% figure;
% showImage(kx, ky, specEx, DisplayFormat="Magnitude");
% 
% figure;
% showImage(kx, ky, specEx_num, DisplayFormat="Magnitude");
% 
% figure;
% showImage(kx, ky, specEx - specEx_num, DisplayFormat="dB");
% 
% figure;
% showImage(kx, ky, specEy, DisplayFormat="Magnitude");
% 
% figure;
% showImage(kx, ky, specEy_num, DisplayFormat="Magnitude");
% 
% figure;
% showImage(kx, ky, specEy - specEy_num, DisplayFormat="dB");

figure;
showImage(x, y, convVal, DisplayFormat="dB", Normalize=true);
clim([-300, 0])

figure;
showImage(x, y, convVal, DisplayFormat="Magnitude", Normalize=true);





