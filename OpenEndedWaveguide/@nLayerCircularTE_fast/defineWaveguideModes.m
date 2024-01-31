function [modeStructs] = defineWaveguideModes(O)
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

%% Define Waveguide ModeStructs
%#ok<*AGROW>
for n = 1:size(modes_TE, 1)
    modeStructs(n, 1) = nLayer.getCircularModeStruct(wgR, 0, n, "TE");
end

end


