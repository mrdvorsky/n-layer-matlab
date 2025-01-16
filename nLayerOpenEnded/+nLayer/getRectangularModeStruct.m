function [modeStruct] = getRectangularModeStruct(m, n, wgA, wgB, TE_TM)
%GETSPECTRUMRECTANGULAR Get function object defining waveguide spectrums.
% This function returns function objects
%

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

%% Define Weighting Functions
if strcmp(TE_TM, "TE")
    signTE = (-1).^ceil(0.5*(m + n)) ...
        .* (-1j).^(m + n - 1);
    scaleTE = signTE .* sqrt(wgA*wgB) ./ (8*pi);
    scaleTE_h = scaleTE * ((n*wgA).^2 + (m*wgB).^2) ./ (m*n*wgA*wgB*kc);
    scaleTE_e = scaleTE * ((n*wgA).^2 - (m*wgB).^2) ./ (m*n*wgA*wgB*kc);

    if m == 0 || n == 0
        scaleTE = scaleTE * sqrt(8);
    end

    if m == 0
        WhSpec = @(kx, ky, kr, kphi) scaleTE * sin(kphi) ...
            .* sinc((0.5/pi).*(wgA.*kx)) .* double_sinc(ky, wgB, n);
        WeSpec = @(kx, ky, kr, kphi) scaleTE * cos(kphi) ...
            .* sinc((0.5/pi).*(wgA.*kx)) .* double_sinc(ky, wgB, n);
    elseif n == 0
        WhSpec = @(kx, ky, kr, kphi) scaleTE * cos(kphi) ...
            .* double_sinc(kx, wgA, m) .* sinc((0.5/pi).*(wgB.*ky));
        WeSpec = @(kx, ky, kr, kphi) -scaleTE * sin(kphi) ...
            .* double_sinc(kx, wgA, m) .* sinc((0.5/pi).*(wgB.*ky));
    else
        WhSpec = @(kx, ky, kr, kphi) scaleTE_h * kr .* sin(2*kphi) ...
            .* double_sinc(kx, wgA, m) .* double_sinc(ky, wgB, n);
        WeSpec = @(kx, ky, kr, kphi) kr .* (scaleTE_e + scaleTE_h.*cos(2*kphi)) ...
            .* double_sinc(kx, wgA, m) .* double_sinc(ky, wgB, n);
    end
else
    signTM = (-1).^ceil(0.5*(m - n + 1)) ...
        .* (-1j).^(m + n + 1);
    scaleTM = signTM .* sqrt(wgA*wgB) ./ (4*pi*kc);

    WhSpec = @(~, ~, ~, ~) 0;
    WeSpec = @(kx, ky, kr, ~) scaleTM * kr .* double_sinc(kx, wgA, m) .* double_sinc(ky, wgB, n);
end

%% Create Spatial Functions
scale_All = 2 ./ ((2*pi) .* kc .* sqrt(wgA .* wgB));
if strcmp(TE_TM, "TE")
    HertzZ = @(x, y) scale_All ...
        .* cos(m.*pi.*(x - 0.5*wgA) ./ wgA) ...
        .* cos(n.*pi.*(y - 0.5*wgB) ./ wgB);
    HertzZ_dx = @(x, y) -(m.*pi./wgA) .* scale_All ...
        .* sin(m.*pi.*(x - 0.5*wgA) ./ wgA) ...
        .* cos(n.*pi.*(y - 0.5*wgB) ./ wgB);
    HertzZ_dy = @(x, y) -(n.*pi./wgB) .* scale_All ...
        .* cos(m.*pi.*(x - 0.5*wgA) ./ wgA) ...
        .* sin(n.*pi.*(y - 0.5*wgB) ./ wgB);
else
    HertzZ = @(x, y) scale_All ...
        .* sin(m.*pi.*(x - 0.5*wgA) ./ wgA) ...
        .* sin(n.*pi.*(y - 0.5*wgB) ./ wgB);
    HertzZ_dx = @(x, y) (m.*pi./wgA) .* scale_All ...
        .* cos(m.*pi.*(x - 0.5*wgA) ./ wgA) ...
        .* sin(n.*pi.*(y - 0.5*wgB) ./ wgB);
    HertzZ_dy = @(x, y) (n.*pi./wgB) .* scale_All ...
        .* sin(m.*pi.*(x - 0.5*wgA) ./ wgA) ...
        .* cos(n.*pi.*(y - 0.5*wgB) ./ wgB);
end

boundaryPoints = [wgA .* [-0.5, 0.5, 0.5, -0.5]; ...
    wgB .* [-0.5, -0.5, 0.5, 0.5]].';

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
    WhSpec=WhSpec, WeSpec=WeSpec, ...
    CutoffWavenumber=kc, ...
    ApertureWidth=hypot(wgA, wgB), ...
    SymmetryX=symmetryX, ...
    SymmetryY=symmetryY, ...
    HertzAz=HertzZ, ...
    HertzAz_dx=HertzZ_dx, ...
    HertzAz_dy=HertzZ_dy, ...
    BoundaryPoints=boundaryPoints);

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

function v = double_sinc(k, a, m)
    v = (                  sinc((0.5/pi) .* (a.*k - m.*pi)) ...
        + (-1).^(m + 1) .* sinc((0.5/pi) .* (a.*k + m.*pi)) );
end




