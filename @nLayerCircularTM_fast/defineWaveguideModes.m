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
wgR = O.waveguideR;

modes_TM = O.modes_TM;

numModes_TE = size(modes_TM, 1);

%% Define Waveguide TE Modes
modeSpecEx_TM = {};
modeSpecEy_TM = {};
cutoffWavenumber_TM = [];
modeScale_TM = [];
modeLabels_TM = strings(0);
%#ok<*AGROW>
for ii = 1:size(modes_TM, 1)
    n = modes_TM(ii, 1);
    modeLabels_TM(ii) = sprintf("TM_{%d,%d}", 0, n);

    [modeSpecEx_TM{ii}, modeSpecEy_TM{ii}, cutoffWavenumber_TM(ii), modeScale_TM(ii)] ...
        = nLayerOpenEnded.getSpectrumCircular(wgR, 0, n, "TM");
end

%% Construct Output Struct
modeStruct = nLayer.createModeStruct(...
    SpecEx_TM=modeSpecEx_TM, ...
    SpecEy_TM=modeSpecEy_TM, ...
    CutoffBeta_TM=cutoffWavenumber_TM, ...
    PhaseScaleFactor_TM=modeScale_TM, ...
    OutputModes_TM=((1:numModes_TE) == 1), ...
    ModeLabels_TM=modeLabels_TM, ...
    ModeSymmetryAxial="TM", ...
    IntegralScaleFactor=2*wgR);

%% Disable Mode Scaling and Orthogonality Check
% The following line can be uncommented after debugging and verifying that
% all modes are orthogonal and are scaled properly, but should be enabled
% during development.

% modeStruct.CheckModeScalingAndOrthogonality = false;

%% Plot Electric Fields
% xPlot = 0.8 * wgA * linspace(-1, 1, 1001);
% yPlot = 0.8 * wgB * linspace(-1, 1, 501);
% nLayerOpenEnded.plotFields(modeStruct, xPlot, yPlot);

end


