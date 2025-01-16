function [modeStruct] = getCoaxialModeStruct(m, n, r_inner, r_outer, TE_TM, isRotated)
%GETCOAXIALMODESTRUCT Get function object defining waveguide spectrums.
% This function returns a modeStruct for the coaxial waveguide modes.

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
% if (m == 0) && (n == 0)
%     scale = 1 ./ sqrt(2*pi * (log(r_outer) - log(r_inner)));
%     Ex = @(kx, ky, kr, kphi) scale .* cos(kphi) .* besseljyTEM(kr, [r_inner, r_outer]);
%     Ey = @(kx, ky, kr, kphi) scale .* sin(kphi) .* besseljyTEM(kr, [r_inner, r_outer]);
% 
%     modeStruct = nLayer.createModeStruct(TE_TM, "TEM", ...
%     ExSpec=Ex, EySpec=Ey, ...
%     WhSpec=@(~, ~, ~, ~) 0, WeSpec=@(~, ~, ~, ~) 0, ...
%     CutoffWavenumber=0, ...
%     ApertureWidth=2*r_outer, ...
%     SymmetryX="PMC", ...
%     SymmetryY="PMC", ...
%     SymmetryAxial="TM");
%     return;
% end

if (m == 0) && (n == 0)
    scaleTEM = 1j ./ sqrt(2*pi * (log(r_outer) - log(r_inner)));

    modeStruct = nLayer.createModeStruct(...
        TE_TM, "TEM", ...
        WhSpec=@(~, ~, ~, ~) 0, ...
        WeSpec=@(~, ~, kr, ~) scaleTEM*(JmOverKr(m, r_outer, kr) - JmOverKr(m, r_inner, kr)), ...
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
else                        % TM
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

%% Define Weighting Functions
if strcmp(TE_TM, "TE")
    signTE = (-1).^((m + 1).*isRotated + m + (m>0) + ceil(0.5*m)) ...
        .* (-1j).^(m + 1);

    J_inner = r_inner*besseljy(a, b, m, r_inner*kc);
    J_outer = r_outer*besseljy(a, b, m, r_outer*kc);

    scaleTE = signTE ./ sqrt(0.5*(1 + (m==0))*pi) ./ kc ./ sqrt(...
        J_outer.^2 .* (1 - (m./(r_outer*kc)).^2) ...
        - J_inner.^2 .* (1 - (m./(r_inner*kc)).^2));

    scaleTE_h_inner = scaleTE .* J_inner .* kc.^2;
    scaleTE_h_outer = scaleTE .* J_outer .* kc.^2;
    scaleTE_e_inner = scaleTE .* J_inner .* m ./ r_inner;
    scaleTE_e_outer = scaleTE .* J_outer .* m ./ r_outer;

    if isRotated
        WhSpec = @(~, ~, kr, kphi) (scaleTE_h_outer*besseljprime(m, r_outer.*kr) - scaleTE_h_inner*besseljprime(m, r_inner.*kr)) ...
            ./ (kr.^2 - kc.^2) .* sin(m*kphi);
        WeSpec = @(~, ~, kr, kphi) -(scaleTE_e_outer*JmOverKr(m, r_outer, kr) - scaleTE_e_inner*JmOverKr(m, r_inner, kr)) ...
           .* cos(m*kphi);
    else
        WhSpec = @(~, ~, kr, kphi) (scaleTE_h_outer*besseljprime(m, r_outer.*kr) - scaleTE_h_inner*besseljprime(m, r_inner.*kr)) ...
            ./ (kr.^2 - kc.^2) .* cos(m*kphi);
        WeSpec = @(~, ~, kr, kphi) (scaleTE_e_outer*JmOverKr(m, r_outer, kr) - scaleTE_e_inner*JmOverKr(m, r_inner, kr)) ...
            .* sin(m*kphi);
    end
    
    if m == 0
        WeSpec = @(~, ~, ~, ~) 0;
    end
else    % TM
    signTM = (-1).^(m*isRotated + 1 + ceil(0.5*m)) ...
        .* (-1j).^(m + 1);

    J_inner = r_inner*besseljyprime(a, b, m, r_inner*kc);
    J_outer = r_outer*besseljyprime(a, b, m, r_outer*kc);

    scaleTM = signTM .* sqrt((2 - (m==0)) ./ pi) ...
        ./ sqrt(J_outer.^2 - J_inner.^2);

    scaleTM_inner = scaleTM .* J_inner;
    scaleTM_outer = scaleTM .* J_outer;

    WhSpec = @(~, ~, ~, ~) 0;
    if isRotated
        WeSpec = @(~, ~, kr, kphi) (scaleTM_outer*besselj(m, r_outer.*kr) - scaleTM_inner*besselj(m, r_inner.*kr)) ...
            .* (kr ./ (kr.^2 - kc.^2)) .* sin(m*kphi);
    else
        WeSpec = @(~, ~, kr, kphi) (scaleTM_outer*besselj(m, r_outer.*kr) - scaleTM_inner*besselj(m, r_inner.*kr)) ...
            .* (kr ./ (kr.^2 - kc.^2)) .* cos(m*kphi);
    end
end

%% Define Symmetries
isX_PEC = true;
isY_PEC = mod(m, 2) == 0;

if xor(strcmp(TE_TM, "TM"), isRotated)  % TM and rotation flip PEC/PMC.
    isX_PEC = ~isX_PEC;
    isY_PEC = ~isY_PEC;
end

% Set symmetry flags.
symmetryX = "PMC";
if isX_PEC
    symmetryX = "PEC";
end

symmetryY = "PMC";
if isY_PEC
    symmetryY = "PEC";
end

symmetryAxial = "None";
if m == 0
    symmetryAxial = TE_TM;
end

%% Create Mode Struct
modeStruct = nLayer.createModeStruct(...
    TE_TM, ...
    sprintf("%s_{%d,%d}", TE_TM, m, n), ...
    ExSpec=ExRot, EySpec=EyRot, ...
    WhSpec=WhSpec, WeSpec=WeSpec, ...
    CutoffWavenumber=kc, ...
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
    scale = (1).^(m) * 0.5 * sqrt(1 + (m~=0)) ...
        ./ (kc .* sqrt(besselj(m, r1r2(2) * kc).^2 - ...
        besselj(m - 1, r1r2(2) * kc).*besselj(m + 1, r1r2(2) * kc)) .* sqrt(pi));
end

function [scale] = scaleFactorTM(r1r2, kc, a, b, m, n)
    int1 = integral(@(x) x .* (besseljy(a, b, m, kc.*x).^2), r1r2(1), r1r2(2));
    
    % scale = (1j).^(m) * 0.5 ./ sqrt(1 + (n~=0)) ...
    %     ./ (kc .* sqrt((besseljyprime(a, b, m, r1r2(2) * kc).^2 ...
    %     - (1).^(n) .* besseljyprime(a, b, m, r1r2(1) * kc).^2)) .* sqrt(pi));

    scale = -(-1).^(m) * 0.5 * sqrt(1 + (m~=0)) ...
        ./ kc .* r1r2(2) ./ sqrt(int1) ./ sqrt(2*pi);
end


function [y] = JmOverKr(m, wgR, kr)
    y = besselj(m, wgR.*kr) ./ kr;
    y(kr == 0) = 0.5 * wgR * (m == 1);
end

