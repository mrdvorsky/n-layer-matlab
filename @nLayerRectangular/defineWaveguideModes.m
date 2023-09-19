function [modeStruct] = defineWaveguideModes(O)
%DEFINEWAVEGUIDEMODES Defines waveguide modes for rectangular waveguide.
% Defines the mode spectrums for a rectangular waveguide. Returns a
% modeStruct as required by the "nLayerOpenEnded" class.
%
% Author: Matt Dvorsky

arguments
    O;
end

%% Get Waveguide Info
wgA = O.waveguideA;
wgB = O.waveguideB;

modes_TE = O.modes_TE;
modes_TM = O.modes_TM;

numModes_TE = size(modes_TE, 1);
numModes_TM = size(modes_TM, 1);

modeSymmetryX = O.modeSymmetryX;
modeSymmetryY = O.modeSymmetryY;

%% Define Waveguide TE Modes
modeSpectrumEx_TE = cell(numModes_TE, 1);
modeSpectrumEy_TE = cell(numModes_TE, 1);
cutoffBeta_TE = zeros(numModes_TE, 1);
for ii = 1:size(modes_TE, 1)
    m = modes_TE(ii, 1);
    n = modes_TE(ii, 2);

    cutoffBeta_TE(ii) = hypot(m*pi/wgA, n*pi/wgB);

    scaleX = -(n/wgB);
    scaleY = (m/wgA);
    scaleBoth = 1 ./ hypot(scaleX, scaleY);

    modeSpectrumEx_TE{ii} = @(kx, ky, ~, ~) scaleX .* scaleBoth ...
        .* C_int(kx, wgA, m) .* S_int(ky, wgB, n);

    modeSpectrumEy_TE{ii} = @(kx, ky, ~, ~) scaleY .* scaleBoth ...
        .* S_int(kx, wgA, m) .* C_int(ky, wgB, n);
end

%% Define Waveguide TM Modes
modeSpectrumEx_TM = cell(numModes_TM, 1);
modeSpectrumEy_TM = cell(numModes_TM, 1);
cutoffBeta_TM = zeros(numModes_TM, 1);
for ii = 1:size(modes_TM, 1)
    m = modes_TM(ii, 1);
    n = modes_TM(ii, 2);

    cutoffBeta_TM(ii) = hypot(m*pi/wgA, n*pi/wgB);

    scaleX = (m/wgA);
    scaleY = (n/wgB);
    scaleBoth = 1 ./ hypot(scaleX, scaleY);

    modeSpectrumEx_TM{ii} = @(kx, ky, ~, ~) scaleX .* scaleBoth ...
        .* C_int(kx, wgA, m) .* S_int(ky, wgB, n);

    modeSpectrumEy_TM{ii} = @(kx, ky, ~, ~) scaleY .* scaleBoth ...
        .* S_int(kx, wgA, m) .* C_int(ky, wgB, n);
end

%% Integral
% int1 = integral(@(k) C_int(k, 1, 0).^2 + 0*k, -inf, inf)
% int1 = integral(@(k) C_int(k, 2, 0).^2 + 0*k, -inf, inf)
% int2 = integral(@(k) C_int(k, 2, 1).^2 + 0*k, -inf, inf)
% int3 = integral(@(k) C_int(k, 1, 2).^2 + 0*k, -inf, inf)
% int4 = integral(@(k) C_int(k, 2, 3).^2 + 0*k, -inf, inf)
% 
% int5 = integral(@(k) S_int(k, 1, 0).^2 + 0*k, -inf, inf)
% int6 = integral(@(k) S_int(k, 2, 1).^2 + 0*k, -inf, inf)
% int7 = integral(@(k) S_int(k, 1, 2).^2 + 0*k, -inf, inf)
% int8 = integral(@(k) S_int(k, 2, 3).^2 + 0*k, -inf, inf)

%% Construct Output Struct
modeStruct = O.createModeStruct(...
    SpecEx_TE=modeSpectrumEx_TE, ...
    SpecEy_TE=modeSpectrumEy_TE, ...
    SpecEx_TM=modeSpectrumEx_TM, ...
    SpecEy_TM=modeSpectrumEy_TM, ...
    CutoffBeta_TE=cutoffBeta_TE, ...
    CutoffBeta_TM=cutoffBeta_TM, ...
    OutputModes_TE=((1:numModes_TE) == 1), ...
    OutputModes_TM=((1:numModes_TM) == 1), ...
    ModeSymmetryX=modeSymmetryX, ...
    ModeSymmetryY=modeSymmetryY, ...
    IntegralScaleFactor=(pi.^2 ./ wgA));

%% Disable Mode Scaling and Orthogonality Check
% The following line can be uncommented after debugging and verifying that
% all modes are orthogonal and are scaled properly, but should be enabled
% during development.

% modeStruct.CheckModeScalingAndOrthogonality = false;

end



%% Integrals over Sin and Cos
function v = S_int(k, a, m)
if m == 0
    v = 0;
    return;
end
v = sqrt(a*pi) * (0.5/pi) ...
    .* (               sinc((0.5/pi) .* (a.*k - m.*pi)) ...
    + (-1).^(m + 1) .* sinc((0.5/pi) .* (a.*k + m.*pi)) );
end


function v = C_int(k, a, m)
if m == 0
    v = sqrt(0.5*a/pi) ...
        .* sinc((0.5/pi) .* (a.*k));
    return;
end
v = a .* sqrt(a/pi) .* (0.5/pi) .* k ./ m ...
    .* (               sinc((0.5/pi) .* (a.*k - m.*pi)) ...
    + (-1).^(m + 1) .* sinc((0.5/pi) .* (a.*k + m.*pi)) );
end

