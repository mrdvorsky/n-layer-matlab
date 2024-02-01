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
wgA = O.waveguideA;
wgB = O.waveguideB;
modes_TE = O.modes_TE;
modes_TM = O.modes_TM;

%% Define Waveguide ModeStructs
modesAll = [modes_TE; modes_TM];
modeTypes = [repmat("TE", size(modes_TE, 1), 1); ...
    repmat("TM", size(modes_TM, 1), 1)];

%#ok<*AGROW>
for ii = 1:size(modesAll, 1)
    m = modesAll(ii, 1);
    n = modesAll(ii, 2);
    modeStructs(ii, 1) = nLayer.getRectangularModeStruct(wgA, wgB, m, n, modeTypes(ii));
end

end


