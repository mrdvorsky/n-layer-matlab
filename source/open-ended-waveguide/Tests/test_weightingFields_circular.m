clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(2, 2, waveguideR=1.2, ...
    modeSymmetryX="None", modeSymmetryY="None");

kx(:, 1) = 10 * linspace(-1, 1, 100);
ky(1, :) = 10 * linspace(-1, 1, 100);

%% Calculate Spectrums
kr = hypot(kx, ky);
kphi = atan2(ky, kx);

for ii = 1:numel(NL.modeStructs)
    Ex = NL.modeStructs(ii).ExSpec;
    Ey = NL.modeStructs(ii).EySpec;
    Wh = NL.modeStructs(ii).WhSpec;
    We = NL.modeStructs(ii).WeSpec;

    specE = cos(kphi).*Ex(kx, ky, kr, kphi) + sin(kphi).*Ey(kx, ky, kr, kphi);
    specH = sin(kphi).*Ex(kx, ky, kr, kphi) - cos(kphi).*Ey(kx, ky, kr, kphi);
    specH2 = Wh(kx, ky, kr, kphi) + 0*kr;
    specE2 = We(kx, ky, kr, kphi);
    
    scaleH(ii) = specH2(:) \ specH(:);
    scaleE(ii) = specE2(:) \ specE(:);

    errH(ii) = max(abs(specH2(:) - specH(:)));
    errE(ii) = max(abs(specE2(:) - specE(:)));

    % figure;
    % showImage(kx, ky, specE);
    % title(NL.modeStructs(ii).ModeLabel);
    % 
    % figure;
    % showImage(kx, ky, specE2);
    % title(NL.modeStructs(ii).ModeLabel);
end

%% Print
for ii = 1:numel(NL.modeStructs)
    fprintf("%s: %.2f, %.2f, %g, %g\n", NL.modeStructs(ii).ModeLabel, ...
        scaleH(ii), scaleE(ii), db(errH(ii)), db(errE(ii)));
end




