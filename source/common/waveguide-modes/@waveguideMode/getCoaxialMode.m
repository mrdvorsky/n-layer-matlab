function [mode] = getCoaxialMode(m, n, r_inner, r_outer, TE_TM, isRotated)
%Get "waveguideMode" object for specific coaxial waveguide mode.
% This function returns a "waveguideMode" for a specific waveguide mode.
%
% For TE modes, the the x-axis will have PEC symmetry, and the y-axis will
% have the same if "m" is even. For TM modes, the the x-axis will have PMC
% symmetry, and the y-axis will have the same if "m" is even. If
% "isRotated" is true, then PMC and PEC symmetries will flip for all
% modes.
%
% Example Usage:
%   [mode] = getCoaxialMode(m, n, wgRi, wgRo, "TE", false);
%   [mode] = getCoaxialMode(m, n, wgRi, wgRo, "TE", true);
%   [mode] = getCoaxialMode(m, n, wgRi, wgRo, "TM", false);
%   [mode] = getCoaxialMode(m, n, wgRi, wgRo, "TM", true);
%
%
% Inputs:
%   m - Scalar value of "m" the returned TEmn or TMmn mode.
%   n - Scalar value of "n" the returned TEmn or TMmn mode.
%   r_inner - Coaxial waveguide inner radius.
%   r_outer - Coaxial waveguide outer radius.
%   TE_TM - String containing "TE" or "TM".
%   isRotated - Whether or not to swap PEC/PMC symmetries, see above.
%
% Author: Matt Dvorsky

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
    scaleTEM = 1j ./ sqrt(2*pi * (log(r_outer) - log(r_inner)));

    mode = waveguideMode(...
        modeLabel="TEM", ...
        modeType=TE_TM, ...
        WhSpec=@(~, ~, ~, ~) 0, ...
        WeSpec=@(~, ~, kr, ~) scaleTEM*(JmOverKr(m, r_outer, kr) - JmOverKr(m, r_inner, kr)), ...
        kc0=0, ...
        apertureSize=2*r_outer, ...
        symmetryX="PMC", ...
        symmetryY="PMC", ...
        symmetryAxial="TM");

    return;
end

%% Calculate Mode Cutoff
if strcmp(TE_TM, "TE")
    [kc, t] = besselCrossPrimeZeros(m, r_outer./r_inner, n);
else    % TM
    [kc, t] = besselCrossZeros(m, r_outer./r_inner, n);
end

kc = kc ./ r_inner;

%% Define Weighting Functions
if strcmp(TE_TM, "TE")
    signTE = (-1).^((m + 1).*isRotated + m + (m>0) + ceil(0.5*m)) ...
        .* (-1j).^(m + 1);

    J_inner = r_inner*besselCylinder(m, t, r_inner*kc);
    J_outer = r_outer*besselCylinder(m, t, r_outer*kc);

    scaleTE = signTE ./ sqrt(0.5*(1 + (m==0))*pi) ./ kc ./ sqrt(...
        J_outer.^2 .* (1 - (m./(r_outer*kc)).^2) ...
        - J_inner.^2 .* (1 - (m./(r_inner*kc)).^2));

    scaleTE_h_inner = scaleTE .* J_inner .* kc.^2;
    scaleTE_h_outer = scaleTE .* J_outer .* kc.^2;
    scaleTE_e_inner = scaleTE .* J_inner .* m ./ r_inner;
    scaleTE_e_outer = scaleTE .* J_outer .* m ./ r_outer;

    if isRotated
        WhSpec = @(~, ~, kr, kphi) (scaleTE_h_outer*besseljPrime(m, r_outer.*kr) ...
            - scaleTE_h_inner*besseljPrime(m, r_inner.*kr)) ...
            ./ (kr.^2 - kc.^2) .* sin(m*kphi);
        WeSpec = @(~, ~, kr, kphi) -(scaleTE_e_outer*JmOverKr(m, r_outer, kr) ...
            - scaleTE_e_inner*JmOverKr(m, r_inner, kr)) ...
            .* cos(m*kphi);
    else
        WhSpec = @(~, ~, kr, kphi) (scaleTE_h_outer*besseljPrime(m, r_outer.*kr) ...
            - scaleTE_h_inner*besseljPrime(m, r_inner.*kr)) ...
            ./ (kr.^2 - kc.^2) .* cos(m*kphi);
        WeSpec = @(~, ~, kr, kphi) (scaleTE_e_outer*JmOverKr(m, r_outer, kr) ...
            - scaleTE_e_inner*JmOverKr(m, r_inner, kr)) ...
            .* sin(m*kphi);
    end

    if m == 0
        WeSpec = @(~, ~, ~, ~) 0;
    end
else    % TM
    signTM = (-1).^(m*isRotated + 1 + ceil(0.5*m)) ...
        .* (-1j).^(m + 1);

    J_inner = r_inner*besselCylinderPrime(m, t, r_inner*kc);
    J_outer = r_outer*besselCylinderPrime(m, t, r_outer*kc);

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
mode = waveguideMode(...
    modeLabel=sprintf("%s_{%d,%d}", TE_TM, m, n), ...
    modeType=TE_TM, ...
    WhSpec=WhSpec, ...
    WeSpec=WeSpec, ...
    kc0=kc, ...
    apertureSize=2*r_outer, ...
    symmetryX=symmetryX, ...
    symmetryY=symmetryY, ...
    symmetryAxial=symmetryAxial);

end


%% Helper Function
function [y] = JmOverKr(m, wgR, kr)
    y = besselj(m, wgR.*kr) ./ kr;
    y(kr == 0) = 0.5 * wgR * (m == 1);
end

