function [modeStruct] = getCircularModeStruct(m, n, wgR, TE_TM, isRotated, options)
%GETCIRCULARMODESTRUCT Get function object defining waveguide spectrums.
% This function returns a modeStruct for the circular waveguide modes.
%
% If TE mode, the the x-axis will be PEC, and the y-axis will be the same
% if "m" is even. For TM modes the the x-axis will be PMC, and the y-axis
% will be the same if "m" is even.
%
% If "isRotated" is true, then PMC and PEC will flip for all modes.
%
% Author: Matt Dvorsky

arguments
    m(1, 1) {mustBeNonnegative, mustBeInteger};
    n(1, 1) {mustBePositive, mustBeInteger};
    wgR(1, 1) {mustBePositive};
    TE_TM(1, 1) string {mustBeMember(TE_TM, ["TE", "TM"])};
    isRotated(1, 1) logical;

    options.kc {mustBeScalarOrEmpty} = [];
end

%% Check Inputs
if isRotated && (m == 0)
    error("Only one TE0n or TM0n mode is possible, but the 'isRotated' " + ...
        "flag was set to true, which will result in a duplicate mode.");
end

%% Get Cutoff Wavenumber
if isempty(options.kc)
    if strcmp(TE_TM, "TE")
        kcAll = besseljprime_zeros(m, n) ./ wgR;
    else
        kcAll = besselj_zeros(m, n) ./ wgR;
    end
    kc = kcAll(end);
else
    kc = options.kc;
end

%% Create Mode Spectrum Functions
if strcmp(TE_TM, "TE")
    scale = scaleFactor_TE(wgR, kc, m);
    if isRotated
        Ex = @(~, ~, kr, kPhi)  circSpectrum_sinSin(kr, kPhi, wgR, kc, m) * scale;
        Ey = @(~, ~, kr, kPhi) -circSpectrum_cosSin(kr, kPhi, wgR, kc, m) * scale;
    else
        Ex = @(~, ~, kr, kPhi)  circSpectrum_sinCos(kr, kPhi, wgR, kc, m) * scale;
        Ey = @(~, ~, kr, kPhi) -circSpectrum_cosCos(kr, kPhi, wgR, kc, m) * scale;
    end
else
    scale = scaleFactor_TM(wgR, kc, m);
    if isRotated
        Ex = @(~, ~, kr, kPhi) -circSpectrum_cosSin(kr, kPhi, wgR, kc, m) * scale;
        Ey = @(~, ~, kr, kPhi) -circSpectrum_sinSin(kr, kPhi, wgR, kc, m) * scale;
    else
        Ex = @(~, ~, kr, kPhi) -circSpectrum_cosCos(kr, kPhi, wgR, kc, m) * scale;
        Ey = @(~, ~, kr, kPhi) -circSpectrum_sinCos(kr, kPhi, wgR, kc, m) * scale;
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
modeStruct = nLayer.createModeStruct(TE_TM, ...
    sprintf("%s_{%d,%d}", TE_TM, m, n), ...
    ExSpec=Ex, EySpec=Ey, ...
    CutoffWavenumber=kc, ...
    ApertureWidth=2*wgR, ...
    SymmetryX=symmetryX, ...
    SymmetryY=symmetryY, ...
    SymmetryAxial=symmetryAxial);

end




%% Helper Functions
function [y] = circSpectrum_cosCos(kr, kphi, wgR, kc, m)
    rot1 = 2 * cos(kphi) .* cos(m .* kphi);
    rot2 = -cos((m + 1) .* kphi);

    int1 = besselj_int1(kr, wgR, kc, m);
    int2 = besselj_int2(kr, wgR, kc, m);

    y = rot1.*int1 + rot2.*int2;
end

function [y] = circSpectrum_sinCos(kr, kphi, wgR, kc, m)
    rot1 = 2 * sin(kphi) .* cos(m .* kphi);
    rot2 = -sin((m + 1) .* kphi);

    int1 = besselj_int1(kr, wgR, kc, m);
    int2 = besselj_int2(kr, wgR, kc, m);

    y = rot1.*int1 + rot2.*int2;
end

function [y] = circSpectrum_cosSin(kr, kphi, wgR, kc, m)
    rot1 = 2 * cos(kphi) .* sin(m .* kphi);
    rot2 = -sin((m + 1) .* kphi);

    int1 = besselj_int1(kr, wgR, kc, m);
    int2 = besselj_int2(kr, wgR, kc, m);

    y = rot1.*int1 + rot2.*int2;
end

function [y] = circSpectrum_sinSin(kr, kphi, wgR, kc, m)
    rot1 = 2 * sin(kphi) .* sin(m .* kphi);
    rot2 = cos((m + 1) .* kphi);

    int1 = besselj_int1(kr, wgR, kc, m);
    int2 = besselj_int2(kr, wgR, kc, m);

    y = rot1.*int1 + rot2.*int2;
end

function [y] = besselj_int1(kr, wgR, kc, m)
    y = kc .* (kr .* besselj(m, wgR.*kr) .* besselj(m - 1, wgR.*kc) ...
        - kc .* besselj(m - 1, wgR.*kr) .* besselj(m, wgR.*kc)) ...
        ./ (kr.^2 - kc.^2);
end

function [y] = besselj_int2(kr, wgR, kc, m)
    JmOverKr = besselj(m, wgR.*kr) ./ kr;
    JmOverKr(kr == 0) = 0.5 * wgR * (m == 1);
    y = (2*m ./ wgR) .* besselj(m, wgR.*kc) .* JmOverKr;
end

function [scale] = scaleFactor_TE(wgR, kc, m)
    scale = 0.5 * sqrt(1 + (m~=0)) ...
        ./ (kc .* sqrt(besselj(m, wgR * kc).^2 - ...
        besselj(m - 1, wgR * kc).*besselj(m + 1, wgR * kc)) .* sqrt(pi));
end

function [scale] = scaleFactor_TM(wgR, kc, m)
    scale = 0.5 * sqrt(1 + (m~=0)) ...
        ./ (kc .* sqrt(besseljprime(m, wgR * kc).^2) .* sqrt(pi));
end



