function [modeStruct] = getRectangularModeStruct(wgA, wgB, m, n, TE_TM)
%GETSPECTRUMRECTANGULAR Get function object defining waveguide spectrums.
% This function returns function objects

arguments
    wgA(1, 1) {mustBePositive};
    wgB(1, 1) {mustBePositive};
    m(1, 1) {mustBeNonnegative, mustBeInteger};
    n(1, 1) {mustBeNonnegative, mustBeInteger};
    TE_TM(1, 1) {mustBeMember(TE_TM, ["TE", "TM"])};
end

%% Create Mode Spectrum Functions
kc = pi * hypot(m/wgA, n/wgB);
scale_All = pi * (1j .^ (m + n + 1)) ./ kc;
if strcmp(TE_TM, "TE")
    scaleX =  (n/wgB) .* scale_All;
    scaleY = -(m/wgA) .* scale_All;
else
    scaleX =  (m/wgA) .* scale_All;
    scaleY =  (n/wgB) .* scale_All;
end

Ex = @(kx, ky, ~, ~) scaleX .* C_int(kx, wgA, m).*S_int(ky, wgB, n);
Ey = @(kx, ky, ~, ~) scaleY .* S_int(kx, wgA, m).*C_int(ky, wgB, n);

%% Create Mode Struct
symmetryX = "Even";
if mod(m, 2) == 0
    symmetryX = "Odd";
end

symmetryY = "Even";
if mod(n, 2) == 0
    symmetryY = "Odd";
end

modeStruct = nLayer.createModeStruct(TE_TM, ...
    sprintf("%s_{%d,%d}", TE_TM, m, n), ...
    ExSpec=Ex, EySpec=Ey, ...
    CutoffWavenumber=kc, MaxOperatingWavenumber=2*kc, ...
    ApertureWidth=hypot(wgA, wgB), ...
    SymmetryX=symmetryX, ...
    SymmetryY=symmetryY);

end




%% Integrals over Sin and Cos
function v = S_int(k, a, m)
if m == 0
    v = 0;
    return;
end
v = sqrt(a*pi) * (0.5/pi) ...
    .* (               sinc((0.5/pi) .* (a.*k - m.*pi)) ...
    + (-1).^(m + 1) .* sinc((0.5/pi) .* (a.*k + m.*pi)) );
end


function v = C_int(k, a, m)
if m == 0
    v = sqrt(0.5*a/pi) ...
        .* sinc((0.5/pi) .* (a.*k));
    return;
end
v = a .* sqrt(a/pi) .* (0.5/pi) .* k ./ m ...
    .* (               sinc((0.5/pi) .* (a.*k - m.*pi)) ...
    + (-1).^(m + 1) .* sinc((0.5/pi) .* (a.*k + m.*pi)) );
end

