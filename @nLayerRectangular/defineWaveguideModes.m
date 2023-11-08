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
%#ok<*AGROW>
for ii = 1:size(modes_TE, 1)
    m = modes_TE(ii, 1);
    n = modes_TE(ii, 2);
    modeLabels_TE(ii) = sprintf("TE_{%d,%d}", m, n);

    [modeSpecEx_TE{ii}, modeSpecEy_TE{ii}, cutoffWavenumber_TE(ii), modeScale_TE(ii)] ...
        = nLayerOpenEnded.getSpectrumRectangular(wgA, wgB, m, n, "TE");
end

%% Define Waveguide TM Modes
for ii = 1:size(modes_TM, 1)
    m = modes_TM(ii, 1);
    n = modes_TM(ii, 2);
    modeLabels_TM(ii) = sprintf("TM_{%d,%d}", m, n);

    [modeSpecEx_TM{ii}, modeSpecEy_TM{ii}, cutoffWavenumber_TM(ii), modeScale_TM(ii)] ...
        = nLayerOpenEnded.getSpectrumRectangular(wgA, wgB, m, n, "TM");
end

%% Construct Output Struct
modeStruct = O.createModeStruct(...
    SpecEx_TE=modeSpecEx_TE, ...
    SpecEy_TE=modeSpecEy_TE, ...
    SpecEx_TM=modeSpecEx_TM, ...
    SpecEy_TM=modeSpecEy_TM, ...
    CutoffBeta_TE=cutoffWavenumber_TE, ...
    CutoffBeta_TM=cutoffWavenumber_TM, ...
    PhaseScaleFactor_TE=modeScale_TE, ...
    PhaseScaleFactor_TM=modeScale_TM, ...
    OutputModes_TE=((1:numModes_TE) == 1), ...
    OutputModes_TM=((1:numModes_TM) == -1), ...
    ModeLabels_TE=modeLabels_TE, ...
    ModeLabels_TM=modeLabels_TM, ...
    ModeSymmetryX=modeSymmetryX, ...
    ModeSymmetryY=modeSymmetryY, ...
    IntegralScaleFactor=(pi.^2 ./ wgA));

%% Disable Mode Scaling and Orthogonality Check
% The following line can be uncommented after debugging and verifying that
% all modes are orthogonal and are scaled properly, but should be enabled
% during development.

% modeStruct.CheckModeScalingAndOrthogonality = false;

%% Plot Electric Fields
xPlot = 0.8 * wgA * linspace(-1, 1, 1001);
yPlot = 0.8 * wgB * linspace(-1, 1, 501);
nLayerOpenEnded.plotFields(modeStruct, xPlot, yPlot);

end


