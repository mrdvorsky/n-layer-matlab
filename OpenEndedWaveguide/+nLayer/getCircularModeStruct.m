function [modeStruct] = getCircularModeStruct(wgR, m, n, TE_TM)
%GETSPECTRUMCircular Get function object defining waveguide spectrums.
% This function returns a modeStruct for the circular waveguide modes.

arguments
    wgR(1, 1) {mustBePositive};
    m(1, 1) {mustBeNonnegative, mustBeInteger};
    n(1, 1) {mustBePositive, mustBeInteger};
    TE_TM(1, 1) {mustBeMember(TE_TM, ["TE", "TM"])};
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
    Ex = @(~, ~, kr, kPhi) -besselIntSin(kr, kPhi, wgR, kc, m) * scale;
    Ey = @(~, ~, kr, kPhi) -besselIntCos(kr, kPhi, wgR, kc, m) * scale;
else
    scale = scaleFactorTM(wgR, kc, m, n);
    Ex = @(~, ~, kr, kPhi) -besselIntCos(kr, kPhi, wgR, kc, m) * scale;
    Ey = @(~, ~, kr, kPhi) +besselIntSin(kr, kPhi, wgR, kc, m) * scale;
end

%% Create Mode Struct
modeStruct = nLayer.createModeStruct(TE_TM, ...
    sprintf("%s_{%d,%d}", TE_TM, m, n), ...
    ExSpec=Ex, EySpec=Ey, ...
    CutoffWavenumber=kc, MaxOperatingWavenumber=2*kc, ...
    ApertureWidth=2*wgR, ...
    SymmetryX="Even", ...
    SymmetryY="Odd");

end




%% Helper Functions
function [y] = besselIntCos(kr, kPhi, wgR, kc, m)
    y = (cos((m - 1) .* kPhi) + cos((m + 1) .* kPhi)) ...
        .* besselInt1(kr, wgR, kc, m) ...
        - cos((m + 1) .* kPhi) .* besselInt2(kr, wgR, kc, m);
end

function [y] = besselIntSin(kr, kPhi, wgR, kc, m)
    y = (sin((m - 1) .* kPhi) - sin((m + 1) .* kPhi)) ...
        .* besselInt1(kr, wgR, kc, m) ...
        + sin((m + 1) .* kPhi) .* besselInt2(kr, wgR, kc, m);
end

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
    scale = (1j).^(m) * 0.5 ./ (kc .* besselj(m, wgR * kc) .* sqrt(pi));
end

function [scale] = scaleFactorTM(wgR, kc, m, n)
    scale = (1j).^(m) * 0.5 ./ (kc .* besseljprime(m, wgR * kc) .* sqrt(pi));
end



