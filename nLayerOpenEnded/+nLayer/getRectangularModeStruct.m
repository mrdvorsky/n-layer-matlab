function [modeStruct] = getRectangularModeStruct(m, n, wgA, wgB, TE_TM)
%GETRECTANGULARMODESTRUCT Get function object defining waveguide spectrums.
% This function returns a modeStruct for the rectangular waveguide modes.
%
% If "m" is even, the y-axis will be PEC, and PMC if odd.
% For "n", the same is true but for the x-axis.
%
% Author: Matt Dvorsky

arguments
    m(1, 1) {mustBeNonnegative, mustBeInteger};
    n(1, 1) {mustBeNonnegative, mustBeInteger};
    wgA(1, 1) {mustBePositive};
    wgB(1, 1) {mustBePositive};
    TE_TM(1, 1) string {mustBeMember(TE_TM, ["TE", "TM"])};
end

%% Check Inputs
if (m == 0) && (n == 0) && strcmp(TE_TM, "TE")
    error("Rectangular mode TE00 does not exist.");
end

if ((m == 0) || (n == 0)) && strcmp(TE_TM, "TM")
    error("Rectangular modes TM0n and TMm0 do not exist.");
end

%% Calculate Mode Cutoff
kc = pi * hypot(m/wgA, n/wgB);

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

%% Create Mode Struct
symmetryY = "PMC";
if mod(m, 2) == 0
    symmetryY = "PEC";
end

symmetryX = "PMC";
if mod(n, 2) == 0
    symmetryX = "PEC";
end

modeStruct = nLayer.waveguideMode(...
    modeLabel=sprintf("%s_{%d,%d}", TE_TM, m, n), ...
    modeType=TE_TM, ...
    WhSpec=WhSpec, ...
    WeSpec=WeSpec, ...
    kc0=kc, ...
    apertureSize=hypot(wgA, wgB), ...
    symmetryX=symmetryX, ...
    symmetryY=symmetryY);

end


%% Helper Function
function v = double_sinc(k, a, m)
    v = (                  sinc((0.5/pi) .* (a.*k - m.*pi)) ...
        + (-1).^(m + 1) .* sinc((0.5/pi) .* (a.*k + m.*pi)) );
end

