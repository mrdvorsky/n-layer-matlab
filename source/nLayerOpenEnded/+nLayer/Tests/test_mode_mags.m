clc;
clear;
close all;

%% Inputs
NL = nLayerRectangular(3, 2, ...
    modeSymmetryX="None", ...
    modeSymmetryY="None", ...
    modeSymmetryAxial="None", ...
    waveguideBand="Ka");

%% Calculate
[kr(:, 1), krw(:, 1)] = fejer2_halfOpen(20000, 2);
[kphi(1, :), kphiw(1, :)] = trap(256, 0, 2*pi);
% kphiw = 4*kphiw;
kx = kr .* cos(kphi);
ky = kr .* sin(kphi);

modeStruct1 = NL.modeStructs(5);
modeStruct2 = NL.modeStructs(9);

We1 = modeStruct1.WeSpec(kx, ky, kr, kphi);
Wh1 = modeStruct1.WhSpec(kx, ky, kr, kphi);
We2 = modeStruct2.WeSpec(-kx, -ky, kr, kphi + pi);
Wh2 = modeStruct2.WhSpec(-kx, -ky, kr, kphi + pi);

Ie = sum(kr .* (We1.*We2) .* kphiw .* krw, "all")
Ih = sum(kr .* (Wh1.*Wh2) .* kphiw .* krw, "all")

Ie + Ih








