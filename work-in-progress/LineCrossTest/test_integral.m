clc;
clear;
close all;

%% Inputs
NL = nLayerRectangular(1, 2, waveguideBand="Ka");
% NL = nLayerCircular(0, 1, waveguideBand="Ka_TE01", modeSymmetryAxial="TM");

%% Integral over Phi
[kphi(1, 1, :), kphi_w(1, 1, :)] = fejer2(800, 0, 0.5*pi);
kphi_w = kphi_w .* (4);
modeFun = @(kr) sum(kphi_w .* ...
    (cos(kphi).*NL.modeStructs(end).ExSpec(kr.*cos(kphi), kr.*sin(kphi), 0*kr, 0*kphi) ...
    + sin(kphi).*NL.modeStructs(end).EySpec(kr.*cos(kphi), kr.*sin(kphi), 0*kr, 0*kphi)).^2, ...
    3);

%% Integral over Kr
[kr, kr_w] = fejer2_halfOpen(32000, 5);

db(1 - sum(modeFun(kr) .* kr .* kr_w))