function [modeStruct] = getCircularModeStruct(m, n, wgR, TE_TM, isRotated)
%GETCIRCULARMODESTRUCT Get function object defining waveguide spectrums.
% This function returns a modeStruct for the circular waveguide modes.

arguments
    m(1, 1) {mustBeNonnegative, mustBeInteger};
    n(1, 1) {mustBePositive, mustBeInteger};
    wgR(1, 1) {mustBePositive};
    TE_TM(1, 1) string {mustBeMember(TE_TM, ["TE", "TM"])};
    isRotated(1, 1) logical;
end

%% Create Mode Spectrum Functions
if strcmp(TE_TM, "TE")
    kc = besseljprime_zeros(m, n) ./ wgR;
else
    kc = besselj_zeros(m, n) ./ wgR;
end
kc = kc(end);

if strcmp(TE_TM, "TE")
    scale = scaleFactorTE(wgR, kc, m, n);
    Ex = @(~, ~, kr, kPhi)  besselIntSin(kr, kPhi, wgR, kc, m) * scale;
    Ey = @(~, ~, kr, kPhi) -besselIntCos(kr, kPhi, wgR, kc, m) * scale;
else
    scale = scaleFactorTM(wgR, kc, m, n);
    Ex = @(~, ~, kr, kPhi) -besselIntCos(kr, kPhi, wgR, kc, m) * scale;
    Ey = @(~, ~, kr, kPhi) -besselIntSin(kr, kPhi, wgR, kc, m) * scale;
end

if ~isRotated
    rotVal = pi/2 ./ max(m, 1);
    ExRot = @(~, ~, kr, kPhi) cos(rotVal) * Ex(0, 0, kr, kPhi - rotVal) ...
        - sin(rotVal) * Ey(0, 0, kr, kPhi - rotVal);
    EyRot = @(~, ~, kr, kPhi) sin(rotVal) * Ex(0, 0, kr, kPhi - rotVal) ...
        + cos(rotVal) * Ey(0, 0, kr, kPhi - rotVal);
else
    ExRot = Ex;
    EyRot = Ey;
end

%% Define Symmetries
isX_PMC = true;
isY_PMC = false;
if mod(m, 2) == 0
    isY_PMC = ~isY_PMC;
end

if xor(strcmp(TE_TM, "TM"), isRotated)
    isX_PMC = ~isX_PMC;
    isY_PMC = ~isY_PMC;
end

symmetryX = "PEC";
if isX_PMC
    symmetryX = "PMC";
end
symmetryY = "PEC";
if isY_PMC
    symmetryY = "PMC";
end
symmetryAxial = "None";
if m == 0
    symmetryAxial = TE_TM;
end

%% Create Mode Struct
modeStruct = nLayer.createModeStruct(TE_TM, ...
    sprintf("%s_{%d,%d}", TE_TM, m, n), ...
    ExSpec=ExRot, EySpec=EyRot, ...
    CutoffWavenumber=kc, MaxOperatingWavenumber=2*kc, ...
    ApertureWidth=2*wgR, ...
    SymmetryX=symmetryX, ...
    SymmetryY=symmetryY, ...
    SymmetryAxial=symmetryAxial);

end




%% Helper Functions
function [y] = besselIntCos(kr, kPhi, wgR, kc, m)
    y = 2 * cos(kPhi) .* cos(m .* kPhi) ...
        .* besselInt1(kr, wgR, kc, m) ...
        - cos((m + 1) .* kPhi) .* besselInt2(kr, wgR, kc, m);
end

function [y] = besselIntSin(kr, kPhi, wgR, kc, m)
    y = 2 * sin(kPhi) .* cos(m .* kPhi) ...
        .* besselInt1(kr, wgR, kc, m) ...
        - sin((m + 1) .* kPhi) .* besselInt2(kr, wgR, kc, m);
end

% function [y] = besselIntCos(kr, kPhi, wgR, kc, m)
%     y = (cos((m - 1) .* kPhi) + cos((m + 1) .* kPhi)) ...
%         .* besselInt1(kr, wgR, kc, m) ...
%         - cos((m + 1) .* kPhi) .* besselInt2(kr, wgR, kc, m);
% end
% 
% function [y] = besselIntSin(kr, kPhi, wgR, kc, m)
%     y = -(sin((m - 1) .* kPhi) - sin((m + 1) .* kPhi)) ...
%         .* besselInt1(kr, wgR, kc, m) ...
%         - sin((m + 1) .* kPhi) .* besselInt2(kr, wgR, kc, m);
% end

function [y] = besselInt1(kr, wgR, kc, m)
    y = kc .* (kr .* besselj(m, wgR.*kr) .* besselj(m - 1, wgR.*kc) ...
        - kc .* besselj(m - 1, wgR.*kr) .* besselj(m, wgR.*kc)) ...
        ./ (kr.^2 - kc.^2);
end

function [y] = besselInt2(kr, wgR, kc, m)
    JmOverKr = besselj(m, wgR.*kr) ./ kr;
    JmOverKr(kr == 0) = 0.5 * wgR * (m == 1);
    y = (2*m ./ wgR) .* besselj(m, wgR.*kc) .* JmOverKr;
end

function [scale] = scaleFactorTE(wgR, kc, m, n)
    scale = (1j).^(m) * 0.5 * sqrt(1 + (m~=0)) ...
        ./ (kc .* sqrt(besselj(m, wgR * kc).^2 - ...
        besselj(m - 1, wgR * kc).*besselj(m + 1, wgR * kc)) .* sqrt(pi));
end

function [scale] = scaleFactorTM(wgR, kc, m, n)
    scale = (1j).^(m) * 0.5 * sqrt(1 + (m~=0)) ...
        ./ (kc .* sqrt(besseljprime(m, wgR * kc).^2) .* sqrt(pi));
end



