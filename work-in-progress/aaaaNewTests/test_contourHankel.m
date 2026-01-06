clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(0, 1, waveguideR=0.5, modeSymmetryAxial="TM");
f = 300;

N = 100;
gam = 0.577215664901532;

%% Modes
specE1 = NL.waveguideModes.WeSpec;

wgR = NL.waveguideR;
kc = NL.mode_kc0;

%% Contour Integration
theta = trap(N, 0, 2*pi);
x1(:, 1) = wgR .* cos(theta);
y1(:, 1) = wgR .* sin(theta);

% x2(1, :) = wgR .* cos(theta + 0.5*mean(diff(theta)));
% y2(1, :) = wgR .* sin(theta + 0.5*mean(diff(theta)));

x2(1, :) = wgR .* 0;
y2(1, :) = wgR .* 1;

h0_mod = @(x) besselh(0, 1, x + eps) - 2j./pi*log(x + eps);
h0_mod = @(x) (besselh(0, 1, x));

r = vectorize(hypot(x1 - x2, y1 - y2), 3);
intFun2 = @(kr) -2*sum(besselj(0, kr .* r), 3) ./ (kr.^2 - kc.^2).^2 .* kr.^3 ...
    ./ N.^1;
intFun3 = @(kr) -2*sum(h0_mod(kr .* r), 3) ./ (kr.^2 - kc.^2).^2 .* kr.^3 ...
    ./ N.^2;


%% Residue
% % 1. Calculate the terms at the pole kc
% % Note: Use besselh(..., 1) for H^(1)
% H0_val = besselh(0, 1, kc .* r);
% H1_val = besselh(1, 1, kc .* r);
% 
% % 2. Compute the Residue of the Double Pole
% % Formula derived from limit( d/dk [ (k-kc)^2 * Integrand ] )
% Residue = -1 ./ (2 * N.^2) .* sum( 2.*H0_val - (kc .* r) .* H1_val, 3 );
% 
% % 3. Apply the Correction
% % The factor is (j * pi) for a pole on the real axis
% Correction = real(1j * pi * Residue);

%% Integral
intFun1 = @(kr) (2*pi) * specE1(0, 0, kr, 0).^2 .* kr;


int1 = integral(intFun1, 0, 30)
int2 = integral(intFun2, 0, 30)
% int3 = integral(intFun3, 0.01, 30*(1 + 1j))

% int3 + Correction


%% Plot
% kr = linspace(0, 50, 1000);

% figure;
% plot(kr, intFun1(kr), "", LineWidth=1.5);
% hold on;
% plot(kr, intFun2(kr), ":", LineWidth=1.5);
% grid on;






