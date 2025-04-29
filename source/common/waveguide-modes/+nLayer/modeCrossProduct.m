function [crossProd] = modeCrossProduct(mode1, mode2, options)
%Calculate mode cross product.
% This function calculates the integral over the cross product of modes.
%
% Author: Matt Dvorsky

arguments
    mode1 nLayer.waveguideMode;
    mode2 nLayer.waveguideMode;

    options.NumIntegralPointsRho(1, 1) {mustBeInteger, mustBePositive} = 8000;
    options.NumIntegralPointsPhi(1, 1) {mustBeInteger, mustBePositive} = 200;
end

%% Calculate Integral Weights and Nodes
Nrho = options.NumIntegralPointsRho;
Nphi = options.NumIntegralPointsPhi;

[kr(:, 1), kr_w(:, 1)] = fejer2_halfOpen(Nrho, 1);
[kphi(1, :), kphi_w(1, :)] = trap(Nphi, 0, 2*pi);
% kphi_w = 4*kphi_w;

%% Calculate Fields
kx = kr .* cos(kphi);
ky = kr .* sin(kphi);

Wh1 = mode1.WhSpec(kx, ky, kr, kphi);
Wh2 = mode2.WhSpec(-kx, -ky, kr, kphi + pi);
We1 = mode1.WeSpec(kx, ky, kr, kphi);
We2 = mode2.WeSpec(-kx, -ky, kr, kphi + pi);

%% Integral
crossProd = -sum(Wh1.*Wh2 .* (kr .* kr_w .* kphi_w), "all") ...
    - sum(We1.*We2 .* (kr .* kr_w .* kphi_w), "all");

end

