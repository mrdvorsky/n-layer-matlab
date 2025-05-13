function [mode] = getCircularMode(m, n, wgR, TE_TM, isRotated)
%Get "waveguideMode" object for specific circular waveguide mode.
% This function returns a "waveguideMode" for a specific waveguide mode.
%
% For TE modes, the the x-axis will have PEC symmetry, and the y-axis will
% have the same if "m" is even. For TM modes, the the x-axis will have PMC
% symmetry, and the y-axis will have the same if "m" is even. If
% "isRotated" is true, then PMC and PEC symmetries will flip for all
% modes.
%
% Example Usage:
%   [mode] = getCircularMode(m, n, wgR, "TE", false);
%   [mode] = getCircularMode(m, n, wgR, "TE", true);
%   [mode] = getCircularMode(m, n, wgR, "TM", false);
%   [mode] = getCircularMode(m, n, wgR, "TM", true);
%
%
% Inputs:
%   m - Scalar value of "m" the returned TEmn or TMmn mode.
%   n - Scalar value of "n" the returned TEmn or TMmn mode.
%   wgR - Circular waveguide radius.
%   TE_TM - String containing "TE" or "TM".
%   isRotated - Whether or not to swap PEC/PMC symmetries, see above.
%
% Author: Matt Dvorsky

arguments
    m(1, 1) {mustBeNonnegative, mustBeInteger};
    n(1, 1) {mustBePositive, mustBeInteger};
    wgR(1, 1) {mustBePositive};
    TE_TM(1, 1) string {mustBeMember(TE_TM, ["TE", "TM"])};
    isRotated(1, 1) logical;
end

%% Check Inputs
if isRotated && (m == 0)
    error("Only one TE0n or TM0n mode is possible, but the 'isRotated' " + ...
        "flag was set to true, which will result in a duplicate mode.");
end

%% Calculate Mode Cutoff
if strcmp(TE_TM, "TE")
    kc = besseljPrimeZeros(m, n) ./ wgR;
else
    kc = besseljZeros(m, n) ./ wgR;
end

%% Define Weighting Functions
if strcmp(TE_TM, "TE")
    signTE = (-1).^((m + 1).*isRotated + m + n + ceil(0.5*m)) ...
        .* (-1j).^(m + 1);
    scaleTE = signTE ./ sqrt(0.5*(1 + (m==0))*pi) ./ kc ./ sqrt(1 - (m./(kc*wgR)).^2);
    scaleTE_h = scaleTE .* kc.^2;
    scaleTE_e = scaleTE .* m ./ wgR;

    if isRotated
        WhSpec = @(~, ~, kr, kphi) scaleTE_h * besseljPrime(m, wgR.*kr) ./ (kr.^2 - kc.^2) .* sin(m*kphi);
        WeSpec = @(~, ~, kr, kphi) -scaleTE_e * JmOverKr(m, wgR, kr) .* cos(m*kphi);
    else
        WhSpec = @(~, ~, kr, kphi) scaleTE_h * besseljPrime(m, wgR.*kr) ./ (kr.^2 - kc.^2) .* cos(m*kphi);
        WeSpec = @(~, ~, kr, kphi) scaleTE_e * JmOverKr(m, wgR, kr) .* sin(m*kphi);
    end

    if m == 0
        WeSpec = @(~, ~, ~, ~) 0;
    end
else    % TM
    signTM = (-1).^(m*isRotated + (m==0) + n + 1 + ceil(0.5*m)) ...
        .* (-1j).^(m + 1);
    scaleTM = signTM .* sqrt((2 - (m==0))/pi);

    WhSpec = @(~, ~, ~, ~) 0;
    if isRotated
        WeSpec = @(~, ~, kr, kphi) scaleTM * besselj(m, wgR.*kr) .* (kr ./ (kr.^2 - kc.^2)) .* sin(m*kphi);
    else
        WeSpec = @(~, ~, kr, kphi) scaleTM * besselj(m, wgR.*kr) .* (kr ./ (kr.^2 - kc.^2)) .* cos(m*kphi);
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
    apertureSize=2*wgR, ...
    symmetryX=symmetryX, ...
    symmetryY=symmetryY, ...
    symmetryAxial=symmetryAxial);

end


%% Helper Function
function [y] = JmOverKr(m, wgR, kr)
    y = besselj(m, wgR.*kr) ./ kr;
    y(kr == 0) = 0.5 * wgR * (m == 1);
end



