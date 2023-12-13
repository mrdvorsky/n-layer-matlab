function [modeStruct] = defineWaveguideModes(O)
%DEFINEWAVEGUIDEMODES Defines waveguide modes for circular waveguide.
% Defines the mode spectrums for a rectangular waveguide. Returns a
% modeStruct as required by the "nLayerOpenEnded" class.
%
% Author: Matt Dvorsky

arguments
    O;
end

%% Get Waveguide Info
wgR = O.waveguideR;

modes_TE = O.modes_TE;
modes_TM = O.modes_TM;

numModes_TE = size(modes_TE, 1);
numModes_TM = size(modes_TM, 1);

modeSymmetryX = O.modeSymmetryX;
modeSymmetryY = O.modeSymmetryY;

%% Define Waveguide TE Modes
%#ok<*AGROW>
for ii = 1:size(modes_TE, 1)
    m = modes_TE(ii, 1);
    n = modes_TE(ii, 2);
    modeLabels_TE(ii) = sprintf("TE_{%d,%d}", m, n);

    [modeSpectrumEx_TE{ii}, modeSpectrumEy_TE{ii}, cutoffBeta_TE(ii)] ...
        = nLayerOpenEnded.getCircularSpectrums(wgR, m, n, "TE");
end

%% Define Waveguide TM Modes
for ii = 1:size(modes_TM, 1)
    m = modes_TM(ii, 1);
    n = modes_TM(ii, 2);
    modeLabels_TM(ii) = sprintf("TM_{%d,%d}", m, n);

    [modeSpectrumEx_TM{ii}, modeSpectrumEy_TM{ii}, cutoffBeta_TM(ii)] ...
        = nLayerOpenEnded.getCircularSpectrums(wgR, m, n, "TM");
end

%% Construct Output Struct
modeStruct = O.createModeStruct(...
    SpecEx_TE=modeSpectrumEx_TE, ...
    SpecEy_TE=modeSpectrumEy_TE, ...
    SpecEx_TM=modeSpectrumEx_TM, ...
    SpecEy_TM=modeSpectrumEy_TM, ...
    CutoffBeta_TE=cutoffBeta_TE, ...
    CutoffBeta_TM=cutoffBeta_TM, ...
    OutputModes_TE=((1:numModes_TE) == 1), ...
    OutputModes_TM=((1:numModes_TM) == -1), ...
    ModeLabels_TE=modeLabels_TE, ...
    ModeLabels_TM=modeLabels_TM, ...
    ModeSymmetryX=modeSymmetryX, ...
    ModeSymmetryY=modeSymmetryY, ...
    IntegralScaleFactor=(pi.^2 ./ wgR));

%% Disable Mode Scaling and Orthogonality Check
% The following line can be uncommented after debugging and verifying that
% all modes are orthogonal and are scaled properly, but should be enabled
% during development.

% modeStruct.CheckModeScalingAndOrthogonality = false;

%% Plot Electric Fields
xPlot = 1.1 * wgR * linspace(-1, 1, 1001);
yPlot = 1.1 * wgR * linspace(-1, 1, 501);
nLayerOpenEnded.plotFields(modeStruct, xPlot, yPlot);

end


