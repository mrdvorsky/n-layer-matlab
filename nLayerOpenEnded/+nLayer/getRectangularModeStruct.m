function [modeStruct] = getRectangularModeStruct(m, n, wgA, wgB, TE_TM)
%GETSPECTRUMRECTANGULAR Get function object defining waveguide spectrums.
% This function returns function objects

arguments
    m(1, 1) {mustBeNonnegative, mustBeInteger};
    n(1, 1) {mustBeNonnegative, mustBeInteger};
    wgA(1, 1) {mustBePositive};
    wgB(1, 1) {mustBePositive};
    TE_TM(1, 1) string {mustBeMember(TE_TM, ["TE", "TM"])};
end

%% Create Mode Spectrum Functions
kc = pi * hypot(m/wgA, n/wgB);
scale_All = pi ./ kc;
if strcmp(TE_TM, "TE")
    scaleX =  (n/wgB) .* scale_All;
    scaleY = -(m/wgA) .* scale_All;
else
    scaleX =  (m/wgA) .* scale_All;
    scaleY =  (n/wgB) .* scale_All;
end

Ex = @(kx, ky, ~, ~) scaleX .* rectSpectrum_cos(kx, wgA, m).*rectSpectrum_sin(ky, wgB, n);
Ey = @(kx, ky, ~, ~) scaleY .* rectSpectrum_sin(kx, wgA, m).*rectSpectrum_cos(ky, wgB, n);

%% Create Mode Struct
symmetryY = "PMC";
if mod(m, 2) == 0
    symmetryY = "PEC";
end

symmetryX = "PMC";
if mod(n, 2) == 0
    symmetryX = "PEC";
end

modeStruct = nLayer.createModeStruct(TE_TM, ...
    sprintf("%s_{%d,%d}", TE_TM, m, n), ...
    ExSpec=Ex, EySpec=Ey, ...
    CutoffWavenumber=kc, ...
    ApertureWidth=hypot(wgA, wgB), ...
    SymmetryX=symmetryX, ...
    SymmetryY=symmetryY);

end




%% Integrals over Sin and Cos
function v = rectSpectrum_sin(k, a, m)
    if m == 0
        v = 0;
        return;
    end
    v = sqrt(a*pi) * (0.5/pi) ...
        .* (               sinc((0.5/pi) .* (a.*k - m.*pi)) ...
        + (-1).^(m + 1) .* sinc((0.5/pi) .* (a.*k + m.*pi)) );
end

function v = rectSpectrum_cos(k, a, m)
    if m == 0
        v = sqrt(0.5*a/pi) ...
            .* sinc((0.5/pi) .* (a.*k));
        return;
    end

    v = a .* sqrt(a/pi) .* (0.5/pi) .* k ./ m ...
        .* (               sinc((0.5/pi) .* (a.*k - m.*pi)) ...
        + (-1).^(m + 1) .* sinc((0.5/pi) .* (a.*k + m.*pi)) );
end




