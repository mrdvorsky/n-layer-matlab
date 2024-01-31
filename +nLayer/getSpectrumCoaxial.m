function [Ex, Ey, cutoffWavenumber, phaseScale] = getSpectrumCoaxial(wgR_inner, wgR_outer, m, n, TE_TM)
%GETSPECTRUMCOAXIAL Get function object defining waveguide spectrums.
% This function returns function objects

arguments
    wgR_inner(1, 1) {mustBePositive};
    wgR_outer(1, 1) {mustBeGreaterThan(wgR_outer, wgR_inner)};
    m(1, 1) {mustBeNonnegative, mustBeInteger};
    n(1, 1) {mustBeNonnegative, mustBeInteger};
    TE_TM(1, 1) {mustBeMember(TE_TM, ["TE", "TM"])};
end

%% Create Mode Spectrum Functions
if (n == 0) && (m == 0)
    if strcmp(TE_TM, "TE")
        error("Coaxial mode TE00 not supported. Use TM00 instead.");
    end
    cutoffWavenumber = 0;
    phaseScale = -1j;

    scale = 1./sqrt(2 * pi * integral(@(kr) besselTEM(kr, wgR_inner, wgR_outer), 0, inf));
    Ex = @(kx, ky, kr, kPhi) scale .* cos(kPhi) .* besselTEM(kr, wgR_inner, wgR_outer);
    Ey = @(kx, ky, kr, kPhi) scale .* sin(kPhi) .* besselTEM(kr, wgR_inner, wgR_outer);
    return;
end

if strcmp(TE_TM, "TE")
    [kc, a, b] = coaxial_cutoffs_TE(wgR_inner, wgR_outer, m, n);
else
    [kc, a, b] = coaxial_cutoffs_TM(wgR_inner, wgR_outer, m, n);
end

% if strcmp(TE_TM, "TE")
%     Ex = @(~, ~, kr, kPhi) -besselIntSin(kr, kPhi, wgR_inner, wgR_outer, kc, a, b, m);
%     Ey = @(~, ~, kr, kPhi) -besselIntCos(kr, kPhi, wgR_inner, wgR_outer, kc, a, b, m);
% else
%     Ex = @(~, ~, kr, kPhi) -besselIntCos(kr, kPhi, wgR_inner, wgR_outer, kc, a, b, m);
%     Ey = @(~, ~, kr, kPhi) besselIntSin(kr, kPhi, wgR_inner, wgR_outer, kc, a, b, m);
% end

if strcmp(TE_TM, "TE")
    Ex = @(~, ~, kr, kPhi) -besselIntSin(kr, kPhi, wgR_inner, wgR_outer, kc, a, b, m);
    Ey = @(~, ~, kr, kPhi) -besselIntCos(kr, kPhi, wgR_inner, wgR_outer, kc, a, b, m);
else
    scale = 1./sqrt(2 * pi * integral(@(kr) kr.^3 .* besselInt(kr, 0, wgR_inner, wgR_outer, kc, a, b, m).^2, 0, inf));
    Ex = @(kx, ky, kr, kPhi) scale .* kx .* besselInt(kr, kPhi, wgR_inner, wgR_outer, kc, a, b, m);
    Ey = @(kx, ky, kr, kPhi) scale .* ky .* besselInt(kr, kPhi, wgR_inner, wgR_outer, kc, a, b, m);
end

%% Set Mode Cutoffs
cutoffWavenumber = kc;

%% Set Phase Scaling Coefficient
phaseScale = (-1j).^(m + 1);

end





%% Temp Bessel Integrals
function [y] = besselTEM(kr, wgR1, wgR2)
    y = (besselj(0, kr .* wgR1) - besselj(0, kr .* wgR2)) ./ kr;
    y(kr == 0) = 0;
end

function [y] = besselInt(kr, kPhi, wgR1, wgR2, kc, a, b, m)
    y = cos(m * kPhi) .* (besselInt1(kr, wgR2, kc, a, b, m) - besselInt1(kr, wgR1, kc, a, b, m));
end

function scale = besselScaleFactor(wgR1, wgR2, kc, a, b, m)
    scale1 = 0.5 * wgR2.^2 .* (besseljy_derivative(a, b, m, kc .* wgR2).^2 ...
        - besseljy_derivative(a, b, m - 1, kc .* wgR2) .* besseljy_derivative(a, b, m + 1, kc .* wgR2));
    scale2 = 0.5 * wgR1.^2 .* (besseljy_derivative(a, b, m, kc .* wgR1).^2 ...
        - besseljy_derivative(a, b, m - 1, kc .* wgR1) .* besseljy_derivative(a, b, m + 1, kc .* wgR1));
    scale = scale1 - scale2;
end

%% Bessel Function Integrals
function [y] = besselIntCos(kr, kPhi, wgR1, wgR2, kc, a, b, m)
    y = (cos((m - 1) .* kPhi) + cos((m + 1) .* kPhi)) ...
        .* (besselInt1(kr, wgR2, kc, a, b, m) - besselInt1(kr, wgR1, kc, a, b, m)) ...
        - cos((m + 1) .* kPhi) ...
        .* (besselInt2(kr, wgR2, kc, a, b, m) - besselInt2(kr, wgR1, kc, a, b, m));
end

function [y] = besselIntSin(kr, kPhi, wgR1, wgR2, kc, a, b, m)
    y = (sin((m - 1) .* kPhi) - sin((m + 1) .* kPhi)) ...
        .* (besselInt1(kr, wgR2, kc, a, b, m) - besselInt1(kr, wgR1, kc, a, b, m)) ...
        + sin((m + 1) .* kPhi) ...
        .* (besselInt2(kr, wgR2, kc, a, b, m) - besselInt2(kr, wgR1, kc, a, b, m));
end

function [y] = besselInt1(kr, wgR, kc, a, b, m)
    y = wgR * (kc .* besselj(m, wgR.*kr) .* besseljy(a, b, m - 1, wgR.*kc) ...
        - kr .* besselj(m - 1, wgR.*kr) .* besseljy(a, b, m, wgR.*kc)) ...
        ./ (kr.^2 - kc.^2);
end

function [y] = besselInt2(kr, wgR, kc, a, b, m)
    JmOverKr = besselj(m, wgR.*kr) ./ kr;
    JmOverKr(kr == 0) = 0.5 * wgR * (m == 1);
    y = (2*m ./ wgR) .* besseljy(a, b, m, wgR.*kc) .* JmOverKr;
end

% function [y] = besselInt1(kr, wgR, kc, a, b, m)
%     y = kc .* (kr .* besselj(m, wgR.*kr) .* besselj(m - 1, wgR.*kc) ...
%         - kc .* besselj(m - 1, wgR.*kr) .* besselj(m, wgR.*kc)) ...
%         ./ (kr.^2 - kc.^2);
% end
% 
% function [y] = besselInt2(kr, wgR, kc, a, b, m)
%     JmOverKr = besselj(m, wgR.*kr) ./ kr;
%     JmOverKr(kr == 0) = 0.5 * wgR * (m == 1);
%     y = (2*m ./ wgR) .* besselj(m, wgR.*kc) .* JmOverKr;
% end

%% Bessel Functions
function [y] = besseljy(a, b, v, x)
    y = a.*besselj(v, x) + b.*bessely(v, x);
end

function [y] = besseljy_derivative(a, b, v, x)
    y = 0.5 * (besseljy(a, b, v - 1, x) - besseljy(a, b, v + 1, x));
end

function [kc, a, b] = coaxial_cutoffs_TM(wgR1, wgR2, m, n)
    kc = besseljy_det_zeros(m, n, wgR1 ./ wgR2) ./ wgR2;
    ab = null([besseljy(1, 0, m, kc .* wgR2), ...
        besseljy(0, 1, m, kc .* wgR2)]);
    a = ab(1);
    b = ab(2);
end

function [kc, a, b] = coaxial_cutoffs_TE(wgR1, wgR2, m, n)
    kc = besseljy_det_derivative_zeros(m, n, wgR1 ./ wgR2) ./ wgR2;
    ab = null([besseljy_derivative(1, 0, m, kc .* wgR2), ...
        besseljy_derivative(0, 1, m, kc .* wgR2)]);
    a = ab(1);
    b = ab(2);
end

%% Bessel Function Zeros
function [x] = besseljy_det_zeros(v, n, a)
    fun = @(y) besseljy_det(v, abs(y), a);
    x = abs(fzero(fun, (v + (v == 0))));
    for ii = 2:n
        x_guess = x + (0.5*pi ./ (1 - a))*[0.9, 1.1];
        while sum(sign(fun(x_guess))) ~= 0
            x_guess(2) = x + (x_guess(2) - x)*1.1;
        end
        x = fzero(fun, x_guess);
    end
end

function [x] = besseljy_det_derivative_zeros(v, n, a)
    fun = @(y) besseljy_derivative_det(v, abs(y), a);
    x = abs(fzero(fun, (v + (v == 0))));
    for ii = 2:n
        x_guess = x + (0.5*pi ./ (1 - a))*[0.9, 1.1];
        while sum(sign(fun(x_guess))) ~= 0
            x_guess(2) = x + (x_guess(2) - x)*1.1;
        end
        x = fzero(fun, x_guess);
    end
end

function [y] = besseljy_det(v, x, a)
    y = besselj(v, x).*bessely(v, a.*x) ...
        - besselj(v, a.*x).*bessely(v, x);
end

function [y] = besseljy_derivative_det(v, x, a)
    y = 0.5 * (besseljy_det(v - 1, x, a) - besseljy_det(v + 1, x, a));
end



