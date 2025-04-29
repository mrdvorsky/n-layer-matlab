clc;
clear;
close all;

%% Inputs
NL = nLayerRectangular(1, 2, waveguideBand="Ka");
modeStruct = NL.modeStructs(end);
kc = modeStruct.CutoffWavenumber;

%% Nodes
[nodes(:, 1), weights(:, 1)] = fejer2(N, -0.5, 0.5);

%% Integration Points
xp(1, :) = modeStruct.BoundaryPoints(:, 1);
yp(1, :) = modeStruct.BoundaryPoints(:, 2);

xps = circshift(xp, -1);
yps = circshift(yp, -1);

xint = (0.5*(xps + xp) + nodes.*(xps - xp));
yint = (0.5*(yps + yp) + nodes.*(yps - yp));
angInt = atan2(yps - yp, xps - xp);

weightsAll = weights .* hypot(xps - xp, yps - yp);

%% Calculate Points
AzN = sin(angInt).*modeStruct.HertzAz_dx(xint, yint) ...
    - cos(angInt).*modeStruct.HertzAz_dy(xint, yint);

%% Calc Spectrum 1
kr = hypot(kx, ky);

specAn1 = -1j*(kx.*modeStruct.ExSpec(kx, ky, 0*kx, 0*ky) ...
    + ky.*modeStruct.EySpec(kx, ky, 0*kx, 0*ky)) ./ kr;

%% Calc Spectrum 2
specAn2 = reshape(...
    nufftn(weightsAll(:).*AzN(:), [xint(:), yint(:)] ./ (2*pi), {kx, ky}), ...
    numel(kx), numel(ky));
specAn2 = specAn2 .* kr ./ (kr.^2 - modeStruct.CutoffWavenumber.^2);

rAll = hypot(xint(:) - xint(:).', yint(:) - yint(:).');
wAll = weightsAll(:).*AzN(:);
wAll = wAll .* wAll.';

[kr_int(1, 1, :), kr_w(1, 1, :)] = fejer2(1000, 0, 30);
% 2*pi*sum(kr_w .* wAll(:) .* besselj(0, kr_int.*rAll(:)) .* kr_int.^3 ./ (kr_int.^2 - kc.^2).^2, "all")

%% Plotting
% figure;
% plot3(real(xint), real(yint), imag(AzN), ".-", LineWidth=1.5);

% figure;
% showImage(kx, ky, specAn1, DisplayFormat="Magnitude");

% figure;
% showImage(kx, ky, specAn2, DisplayFormat="Magnitude");


% figure;
% showImage(kx, ky, (specAn1 - specAn2) ./ max(abs(specAn1(:)), [], "all"), DisplayFormat="dB");
% clim([-100, 0]);





