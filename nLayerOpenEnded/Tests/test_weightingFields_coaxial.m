clc;
clear;
close all;

%% Inputs
NL = nLayerCoaxial(3, 3, waveguideRi=1.3, waveguideRo=2.1, ...
    modeSymmetryAxial="None", modeSymmetryX="PEC", modeSymmetryY="None");

kx(:, 1) = 10 * linspace(-1, 1, 100);
ky(1, :) = 10 * linspace(-1, 1, 100);

%% Calculate Spectrums
% kr = hypot(kx, ky);
% kphi = atan2(ky, kx);

% for ii = 1:numel(NL.modeStructs)
%     Ex = NL.modeStructs(ii).ExSpec;
%     Ey = NL.modeStructs(ii).EySpec;
%     Wh = NL.modeStructs(ii).WhSpec;
%     We = NL.modeStructs(ii).WeSpec;
% 
%     specE = cos(kphi).*Ex(kx, ky, kr, kphi) + sin(kphi).*Ey(kx, ky, kr, kphi);
%     specH = sin(kphi).*Ex(kx, ky, kr, kphi) - cos(kphi).*Ey(kx, ky, kr, kphi);
%     specH2 = Wh(kx, ky, kr, kphi) + 0*kr;
%     specE2 = We(kx, ky, kr, kphi) + 0*kr;
% 
%     scaleH(ii) = specH2(:) \ specH(:);
%     scaleE(ii) = specE2(:) \ specE(:);
% 
%     errH(ii) = max(abs(specH2(:) - specH(:)));
%     errE(ii) = max(abs(specE2(:) - specE(:)));
% 
%     % if strcmp(NL.modeStructs(ii).ModeType, "TM")
%     %     % specE = specE .* (kr.^2 - NL.modeStructs(ii).CutoffWavenumber.^2) ./ kr;
%     %     figure;
%     %     showImage(kx, ky, specE);
%     %     title(NL.modeStructs(ii).ModeLabel);
%     % 
%     %     figure;
%     %     showImage(kx, ky, specE2);
%     %     title(NL.modeStructs(ii).ModeLabel);
%     % end
% end

%% Integrals
[int_kr(:, 1), int_krw(:, 1)] = fejer2_halfOpen(80000, 1);
[int_kphi(1, :), int_kphiw(1, :)] = trap(16, 0, 2*pi);
int_kx = int_kr .* cos(int_kphi);
int_ky = int_kr .* sin(int_kphi);

for ii = 1:numel(NL.modeStructs)
    Ex = NL.modeStructs(ii).ExSpec;
    Ey = NL.modeStructs(ii).EySpec;
    Wh = NL.modeStructs(ii).WhSpec;
    We = NL.modeStructs(ii).WeSpec;

    specEx = Ex(int_kx, int_ky, int_kr, int_kphi);
    specEy = Ey(int_kx, int_ky, int_kr, int_kphi);

    specWe = We(int_kx, int_ky, int_kr, int_kphi);
    specWh = Wh(int_kx, int_ky, int_kr, int_kphi);

    % specE_mag = sum((specEx.^2 + specEy.^2) .* int_kr .* int_kphiw .* int_krw, "all");
    specE_mag = sum((specWe.^2) .* int_kr .* int_kphiw .* int_krw, "all");
    specE_mag = sum((specWe.^2 + specWh.^2) .* int_kr .* int_kphiw .* int_krw, "all");

    kc = NL.modeStructs(ii).CutoffWavenumber;
    if strcmp(NL.modeStructs(ii).ModeType, "TM") || 1
        fprintf("%s: %.4f\n", NL.modeStructs(ii).ModeLabel, specE_mag);
    end
end


%% Print
% for ii = 1:numel(NL.modeStructs)
%     if strcmp(NL.modeStructs(ii).ModeType, "TM")
%         fprintf("%s: %.2f, %.2f, %g, %g\n", NL.modeStructs(ii).ModeLabel, ...
%             (scaleH(ii)), (scaleE(ii)), db(errH(ii)), db(errE(ii)));
%     end
% end




