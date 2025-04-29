clc;
clear;
close all;

%% Inputs
NL = nLayerRectangular(1, 2, waveguideBand="Ka");
% NL = nLayerCircular(0, 1, waveguideBand="Ka_TE01", modeSymmetryAxial="TM");
modeStruct = NL.modeStructs(end);

%% Integral over Phi Analytical
[kphi(1, 1, :), kphi_w(1, 1, :)] = fejer2(800, 0, 2*pi);
kphi_w = kphi_w .* (4);
modeFun = @(kr) sum(kphi_w .* ...
    (cos(kphi).*NL.modeStructs(end).ExSpec(kr.*cos(kphi), kr.*sin(kphi), 0*kr, 0*kphi) ...
    + sin(kphi).*NL.modeStructs(end).EySpec(kr.*cos(kphi), kr.*sin(kphi), 0*kr, 0*kphi)).^2, ...
    3);

%% Integral over Kr
% [kr, kr_w] = fejer2_halfOpen(4*16000, 2);
[kr, kr_w] = fejer2(2*16000, 0 + 1j, 1000 + 1j);
[kr2, kr_w2] = fejer2(1000, 1000 + 1j, 1000);
[kr3, kr_w3] = fejer2(1000, 0, 1j);

kr = [kr; kr2; kr3];
kr_w = [kr_w; kr_w2; kr_w3];

err_db = db(1 - sum(modeFun(kr) .* kr .* kr_w))






