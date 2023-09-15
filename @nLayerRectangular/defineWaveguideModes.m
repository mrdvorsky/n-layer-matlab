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

%% Define Waveguide TE Modes
modeSpectrumEx_TE = cell(numModes_TE, 1);
modeSpectrumEy_TE = cell(numModes_TE, 1);
cutoffBeta_TE = zeros(numModes_TE, 1);
for ii = 1:size(modes_TE, 1)
    m = modes_TE(ii, 1);
    n = modes_TE(ii, 2);

    cutoffBeta_TE(ii) = hypot(m*pi/wgA, n*pi/wgB);

    scaleFactor = pi / (cutoffBeta_TE(ii) * sqrt(wgA*wgB));

    modeSpectrumEx_TE{ii} = @(kx, ky, ~, ~) -(n/wgB) * scaleFactor ...
        .* C_int(kx, wgA, m) .* S_int(ky, wgB, n);

    modeSpectrumEy_TE{ii} = @(kx, ky, ~, ~) (m/wgA) * scaleFactor ...
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

    scaleFactor = pi / (cutoffBeta_TM(ii) * sqrt(wgA*wgB));

    modeSpectrumEx_TM{ii} = @(kx, ky, ~, ~) (m/wgA) * scaleFactor ...
        .* C_int(kx, wgA, m) .* S_int(ky, wgB, n);

    modeSpectrumEy_TM{ii} = @(kx, ky, ~, ~) (n/wgB) * scaleFactor ...
        .* S_int(kx, wgA, m) .* C_int(ky, wgB, n);
end

%% Construct Output Struct
modeStruct = O.createModeStruct(...
    SpecEx_TE=modeSpectrumEx_TE, ...
    SpecEy_TE=modeSpectrumEy_TE, ...
    SpecEx_TM=modeSpectrumEx_TM, ...
    SpecEy_TM=modeSpectrumEy_TM, ...
    CutoffBeta_TE=cutoffBeta_TE, ...
    CutoffBeta_TM=cutoffBeta_TM, ...
    ModeSymmetryX="None", ModeSymmetryY="None", ...
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

v = sqrt( pi ) * m .* sinc((0.5/pi) .* (a.*k - pi*m)) ./ (k  + m.*pi./a);
% v = sqrt( pi ) * m .* (0.25./pi) ...
%     .* (sinc((0.5/pi) .* (a.*k - pi*m)) - sinc((0.5/pi) .* (a.*k + pi*m)));

end

function v = C_int(k, a, m)

v = sqrt(1/pi)*a * k .* sinc((0.5/pi) .* (a.*k - pi*m)) ./ (k  + m.*pi./a);
if m == 0
    v = v .* sqrt(0.5);
end

end

