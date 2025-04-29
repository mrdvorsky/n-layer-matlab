clc;
clear;
close all;

%% Inputs
NL = nLayerRectangular(1, 2);
modeStruct = NL.modeStructs(end);
kc = modeStruct.CutoffWavenumber;

ExSpec = modeStruct.ExSpec;
EySpec = modeStruct.EySpec;
W_TM = @(kr, kphi) (cos(kphi).*ExSpec(kr.*cos(kphi), kr.*sin(kphi), 0, 0) ...
    + sin(kphi).*EySpec(kr.*cos(kphi), kr.*sin(kphi), 0, 0)).^2;

%% Integration 1
[kphi(1, 1, :), kphi_w(1, 1, :)] = fejer2(1*128, 0, pi/2);
kphi_w = 4 * kphi_w;

W_kr1 = @(kr) sum(W_TM(kr, kphi) .* kphi_w, 3);

%% Integral over kr
[kr, kr_w] = fejer2_halfOpen(3000, 2*kc);

int1 = integral(W_kr1, 0, inf, RelTol=1e-8)
int2 = sum(W_kr1(kr) .* kr_w)

err_db = db((int1 - int2) ./ int1)













