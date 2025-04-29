clc;
clear;
close all;

%% Inputs
m = 3;
n = 3;

L = 4;

NL = nLayerRectangular(m, n, waveguideBand="ka", modeSymmetryX="None", modeSymmetryY="None");
wgA = NL.waveguideA;
wgB = NL.waveguideB;

mode1 = NL.modeStructs(find(strcmp(NL.modeLabels, "TE_{1,0}")));
mode2 = NL.modeStructs(find(strcmp(NL.modeLabels, "TE_{1,0}")));

[kphi(1, 1, :), kphi_w(1, 1, :)] = fejer2(200, 0, pi/2);
kphi_w = 4 * kphi_w;

%% Analytical
spec1Fun = @(kr, kphi) ...
    cos(kphi) .* mode1.ExSpec(kr .* cos(kphi), kr .* sin(kphi)) ...
    + sin(kphi) .* mode1.EySpec(kr .* cos(kphi), kr .* sin(kphi));
spec2Fun = @(kr, kphi) ...
    cos(kphi) .* mode2.ExSpec(-kr .* cos(kphi), -kr .* sin(kphi)) ...
    + sin(kphi) .* mode2.EySpec(-kr .* cos(kphi), -kr .* sin(kphi));



specIntFun = @(kr) printSize(kr) + kr .* sum(spec1Fun(kr, kphi) .* spec2Fun(kr, kphi), 3) ./ (L + kr);

Lch = 1;
Lcw = 10;
krToKrc = @(x) x + 1j * (Lch.*Lcw) * x ./ sqrt((Lch + x.^2).*(3*Lcw.^2 + x.^2));
krToKrcPrime = @(x) 1 + 1j * (Lch.*Lcw) * (3*Lch*Lcw.^2 - x.^4) ./ ((Lch + x.^2).*(3*Lcw.^2 + x.^2)).^1.5;

int1 = integral(specIntFun, 0, inf, AbsTol=1e-5, RelTol=0);
int2 = integral(@(kr) specIntFun(krToKrc(kr)) .* krToKrcPrime(kr), 0, inf, AbsTol=1e-3, RelTol=0);


%%






%% Helper
function v = printSize(kr)
    fprintf("%d, ", numel(kr));
    v = 0;
end


