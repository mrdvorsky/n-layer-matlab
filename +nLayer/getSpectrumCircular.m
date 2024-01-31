function [Ex, Ey, cutoffWavenumber, phaseScale] = getSpectrumCircular(wgR, m, n, TE_TM)
%GETSPECTRUMCircular Get function object defining waveguide spectrums.
% This function returns function objects

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
    Ex = @(~, ~, kr, kPhi) -besselIntSin(kr, kPhi, wgR, kc, m);
    Ey = @(~, ~, kr, kPhi) -besselIntCos(kr, kPhi, wgR, kc, m);
else
    Ex = @(~, ~, kr, kPhi) -besselIntCos(kr, kPhi, wgR, kc, m);
    Ey = @(~, ~, kr, kPhi) besselIntSin(kr, kPhi, wgR, kc, m);
end

%% Set Mode Cutoffs
cutoffWavenumber = kc;

%% Set Phase Scaling Coefficient
phaseScale = (-1j).^(m + 1);

end




%% Bessel Function Zeros
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



