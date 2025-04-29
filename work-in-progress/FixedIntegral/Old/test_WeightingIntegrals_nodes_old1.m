clc;
clear;
close all;

%% Inputs
Lc = 0.1;
ccL = 0.3;

ccOrder = 10;

N1 = 400;
N2 = 10;

% krComp(:, 1) = linspace(0, 80, 1001);
% compAngDeg = 70;

%% nLayer
NL = nLayerRectangular(1, 0, waveguideBand="x", convergenceAbsTol=1e-7);

%% Get Spectrum Functions
[specFunEx, specFunEy] = nLayerOpenEnded.getSpectrumRectangular(...
    NL.waveguideA, NL.waveguideB, 1, 0, "TE");

%% Helper
phi1(1, 1, 1, :) = 2*pi * ((1:4*N1) - 0.5) ./ (4*N1);
phi1_weights = 0*phi1(1:N1) + 4*(2*pi ./ (4*N1));
phi1 = phi1(1:N1);

% modeFunH_kr = @(kr, ph) sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
%     .* cos(phi1).^2 .* phi1_weights, 4);
% modeFunH = @(kr) kr .* modeFunH_kr(kr);
modeFunE_kr = @(kr) sum(specFunEy(kr .* cos(phi1), kr .* sin(phi1)).^2 ...
    .* sin(phi1).^2 .* phi1_weights, 4);
% modeFunE = @(kr) kr .* modeFunE_kr(kr);

xToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + x);
xToKrc_weights = @(x) 1 + 1j .* Lc.^2 ./ (x + Lc).^2;

%% Helper 2
% phi2(1, :) = 2*pi * ((1:4*N2) - 0.5) ./ (4*N2);
% phi2_weights = 0*phi2 + 1*(2*pi ./ (4*N2));
% % phi2 = phi2(1:N2);

[phi2(1, :), phi2_weights(1, :)] = fejer2(N2, 0, pi/2);
phi2_weights = 4*phi2_weights;

modeFunE_krphi = @(kr, ph) specFunEy(kr .* cos(ph), kr .* sin(ph)).^2 ...
    .* sin(ph).^2;

%% Compute Weighting 1
[x, x_weights_E1] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
    WeightingFunction=@(x) modeFunE_kr(xToKrc(x)) .* xToKrc_weights(x));

%% Compute Weighting 2
for pp = 1:numel(phi2)
    w = @(x) modeFunE_krphi(xToKrc(x), phi2(pp)) .* xToKrc_weights(x);
    % [~, x_weights_Ephi(:, pp)] = clenshawCurtisHalfOpen(ccOrder, ccL, ...
    %     WeightingFunction=w);

    for cc = 0:ccOrder
        moments(cc + 1, 1) = 0.5 ./ ccL .* integral(...
            @(x) w(x) ...
            .* cos(2*cc * acot(sqrt(x ./ ccL))), ...
            0, inf);
    end
    weights = ifft([moments; moments(end - 1:-1:2)]);
    weights = weights(1:ccOrder + 1);
    weights(2:end - 1) = 2*weights(2:end - 1);

    x_weights_Ephi(:, pp) = 2*ccL .* weights(2:end);
end

x_weights_E2 = sum(x_weights_Ephi .* phi2_weights, 2)


errdB2 = db(x_weights_E1 - x_weights_E2)

%% Compute Weighting 3
rSamp(1, 1, :) = linspace(-20, 20, 1001);
krSamp = fftCoordinates(rSamp, ApplyFftShift=false);
dr = abs(rSamp(2) - rSamp(1));

for pp = 1:numel(phi2)
    pp
    convFun(1, 1, :) = fftshift(ifft(modeFunE_krphi(krSamp, phi2(pp)))) ...
        .* numel(krSamp) .* abs(krSamp(2) - krSamp(1));

    modeFunE_krphi3 = @(kr) (1./(2*pi)) .* dr .* sum(convFun .* exp(-1j .* rSamp .* kr), 3);

    w = @(x) modeFunE_krphi3(xToKrc(x)) .* xToKrc_weights(x);

    for cc = 0:ccOrder
        moments(cc + 1, 1) = 0.5 ./ ccL .* integral(...
            @(x) w(x) ...
            .* cos(2*cc * acot(sqrt(x ./ ccL))), ...
            0, 100);
    end
    weights = ifft([moments; moments(end - 1:-1:2)]);
    weights = weights(1:ccOrder + 1);
    weights(2:end - 1) = 2*weights(2:end - 1);

    x_weights_Ephi(:, pp) = 2*ccL .* weights(2:end);
end

x_weights_E3 = sum(x_weights_Ephi .* phi2_weights, 2)


errdB3 = db(x_weights_E1 - x_weights_E3)

%% Test Fourier 1D
convFunAdj = convFun(rSamp > -inf);
rSampAdj = rSamp(rSamp > -inf);

% convFunAdj = convFun(rSamp > 0);
% rSampAdj = rSamp(rSamp > 0);
convFunAdj = convFun(rSamp <= 0);
rSampAdj = rSamp(rSamp <= 0);

convFunAdj = convFun(501);
rSampAdj = rSamp(501);

modeFunE_krphi3 = @(kr) (1./(2*pi)) .* dr .* sum(convFunAdj .* exp(-1j .* rSampAdj .* kr), 3);
w = @(x) modeFunE_krphi3(xToKrc(x)) .* xToKrc_weights(x);

kxPlot(:, 1) = linspace(0, 10, 501);
kyPlot(1, :) = linspace(-1, 1, 501);

Img = w(kxPlot + 1j*kyPlot);

%% Plot
figure;
showImage(kxPlot, kyPlot, Img, DisplayFormat="dB");
axis normal;
clim([-100, 50]);

figure;
showImage(kxPlot, kyPlot, Img, DisplayFormat="Phase");
axis normal;
colormap hsv;
% clim([-100, 50]);

figure;
plot(kxPlot, db(Img(:, kyPlot == 0)), "", LineWidth=1.5);
grid on;

figure;
plot(kyPlot, db(Img(kxPlot == 0, :)), "", LineWidth=1.5);
grid on;

