function [modeStruct] = getCoaxialModeStruct(m, n, r_inner, r_outer, TE_TM, isRotated)
%GETCOAXIALMODESTRUCT Get function object defining waveguide spectrums.
% This function returns a modeStruct for the circular waveguide modes.

arguments
    m(1, 1) {mustBeNonnegative, mustBeInteger};
    n(1, 1) {mustBeNonnegative, mustBeInteger};
    r_inner(1, 1) {mustBePositive};
    r_outer(1, 1) {mustBePositive, mustBeGreaterThan(r_outer, r_inner)};
    TE_TM(1, 1) string {mustBeMember(TE_TM, ["TE", "TM"])};
    isRotated(1, 1) logical;
end

%% Check Inputs
if (m == 0) && (n == 0) && strcmp(TE_TM, "TE")
    error("Coaxial mode TE00 not supported (to avoid duplicate " + ...
        "modes). Use TM00 instead.");
end

if (m == 0) && isRotated
    error("The 'isRotated' parameter should be false for TE0n " + ...
        "and TM0n modes (to avoid duplicate modes).");
end

%% Handle TEM Mode
if (m == 0) && (n == 0)
    scale = 1 ./ sqrt(2*pi * (log(r_outer) - log(r_inner)));
    Ex = @(kx, ky, kr, kphi) scale .* cos(kphi) .* besseljyTEM(kr, [r_inner, r_outer]);
    Ey = @(kx, ky, kr, kphi) scale .* sin(kphi) .* besseljyTEM(kr, [r_inner, r_outer]);

    modeStruct = nLayer.createModeStruct(TE_TM, "TEM", ...
    ExSpec=Ex, EySpec=Ey, ...
    CutoffWavenumber=0, ...
    ApertureWidth=2*r_outer, ...
    SymmetryX="PMC", ...
    SymmetryY="PMC", ...
    SymmetryAxial="TM");
    return;
end

%% Create Mode Spectrum Functions
if strcmp(TE_TM, "TE")
    [kc, alpha, beta] = besseljyprime_zeros(m, n, r_inner, r_outer);
else
    [kc, alpha, beta] = besseljy_zeros(m, n, r_inner, r_outer);
end
kc = kc(end);
a = alpha(end);
b = beta(end);

r1r2 = [r_inner, r_outer];

if strcmp(TE_TM, "TE")
    scale = scaleFactorTE(r1r2, kc, a, b, m, n);
    Ex = @(~, ~, kr, kPhi)  besseljyIntSin(kr, kPhi, r1r2, kc, a, b, m) * scale;
    Ey = @(~, ~, kr, kPhi) -besseljyIntCos(kr, kPhi, r1r2, kc, a, b, m) * scale;
else
    scale = scaleFactorTM(r1r2, kc, a, b, m, n);
    Ex = @(~, ~, kr, kPhi) -besseljyIntCos(kr, kPhi, r1r2, kc, a, b, m) * scale;
    Ey = @(~, ~, kr, kPhi) -besseljyIntSin(kr, kPhi, r1r2, kc, a, b, m) * scale;
end

if isRotated
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
    ApertureWidth=2*r_outer, ...
    SymmetryX=symmetryX, ...
    SymmetryY=symmetryY, ...
    SymmetryAxial=symmetryAxial);

end




%% Helper Functions
function [y] = besseljyTEM(kr, r1r2)
    y = (besselj(0, kr .* r1r2(1)) - besselj(0, kr .* r1r2(2))) ./ kr;
    y(kr == 0) = 0;
end

function [y] = besseljyIntCos(kr, kPhi, r1r2, kc, a, b, m)
    Int1 = besseljyInt1(kr, r1r2(2), kc, a, b, m) ...
        - (r1r2(1)/r1r2(2)) * besseljyInt1(kr, r1r2(1), kc, a, b, m);
    Int2 = besseljyInt2(kr, r1r2(2), kc, a, b, m) ...
        - (r1r2(1)/r1r2(2)) * besseljyInt2(kr, r1r2(1), kc, a, b, m);

    rad1 = 2 * cos(kPhi) .* cos(m .* kPhi);
    rad2 = -cos((m + 1) .* kPhi);
    
    y = Int1.*rad1 + Int2.*rad2;
end

function [y] = besseljyIntSin(kr, kPhi, r1r2, kc, a, b, m)
    Int1 = besseljyInt1(kr, r1r2(2), kc, a, b, m) ...
        - (r1r2(1)/r1r2(2)) * besseljyInt1(kr, r1r2(1), kc, a, b, m);
    Int2 = besseljyInt2(kr, r1r2(2), kc, a, b, m) ...
        - (r1r2(1)/r1r2(2)) * besseljyInt2(kr, r1r2(1), kc, a, b, m);

    rad1 = 2 * sin(kPhi) .* cos(m .* kPhi);
    rad2 = -sin((m + 1) .* kPhi);
    
    y = Int1.*rad1 + Int2.*rad2;
end

function [y] = besseljyInt1(kr, r, kc, a, b, m)
    y = kc .* (kr .* besselj(m, r.*kr) .* besseljy(a, b, m - 1, r.*kc) ...
        - kc .* besselj(m - 1, r.*kr) .* besseljy(a, b, m, r.*kc)) ...
        ./ (kr.^2 - kc.^2);
end

function [y] = besseljyInt2(kr, r, kc, a, b, m)
    JmOverKr = besselj(m, r.*kr) ./ kr;
    JmOverKr(kr == 0) = 0.5 * r * (m == 1);
    y = (2*m ./ r) .* besseljy(a, b, m, r.*kc) .* JmOverKr;
end

function [scale] = scaleFactorTE(r1r2, kc, a, b, m, n)
    scale = (1j).^(m) * 0.5 * sqrt(1 + (m~=0)) ...
        ./ (kc .* sqrt(besselj(m, r1r2(2) * kc).^2 - ...
        besselj(m - 1, r1r2(2) * kc).*besselj(m + 1, r1r2(2) * kc)) .* sqrt(pi));
end

function [scale] = scaleFactorTM(r1r2, kc, a, b, m, n)
    int1 = integral(@(x) x .* (besseljy(a, b, m, kc.*x).^2), r1r2(1), r1r2(2));
    
    % scale = (1j).^(m) * 0.5 ./ sqrt(1 + (n~=0)) ...
    %     ./ (kc .* sqrt((besseljyprime(a, b, m, r1r2(2) * kc).^2 ...
    %     - (1).^(n) .* besseljyprime(a, b, m, r1r2(1) * kc).^2)) .* sqrt(pi));

    scale = (1j).^(m) * 0.5 ...
        ./ kc .* r1r2(2) ./ sqrt(int1) ./ sqrt(2*pi);
end



