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

modes_TE = O.modes_TE;

numModes_TE = size(modes_TE, 1);

%% Define Waveguide TE Modes
modeSpecEx_TE = {};
modeSpecEy_TE = {};
cutoffWavenumber_TE = [];
modeScale_TE = [];
modeLabels_TE = strings(0);
%#ok<*AGROW>
for ii = 1:size(modes_TE, 1)
    n = modes_TE(ii, 1);
    modeLabels_TE(ii) = sprintf("TE_{%d,%d}", 0, n);

    [modeSpecEx_TE{ii}, modeSpecEy_TE{ii}, cutoffWavenumber_TE(ii), modeScale_TE(ii)] ...
        = nLayerOpenEnded.getSpectrumCircular(wgR, 0, n, "TE");
end

%% Construct Output Struct
modeStruct = nLayer.createModeStruct(...
    SpecEx_TE=modeSpecEx_TE, ...
    SpecEy_TE=modeSpecEy_TE, ...
    CutoffBeta_TE=cutoffWavenumber_TE, ...
    PhaseScaleFactor_TE=modeScale_TE, ...
    OutputModes_TE=((1:numModes_TE) == 1), ...
    ModeLabels_TE=modeLabels_TE, ...
    ModeSymmetryAxial="TE", ...
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


