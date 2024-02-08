function [modeStructs] = defineWaveguideModes(O, symmetryX, symmetryY, symmetryAxial)
%DEFINEWAVEGUIDEMODES Defines waveguide modes for circular waveguide.
% Defines the mode spectrums for a circular waveguide. Returns a
% modeStruct as required by the "nLayerOpenEnded" class.
%
% Author: Matt Dvorsky

arguments
    O;
    symmetryX string {mustBeMember(symmetryX, ["PEC", "PMC", "None"])};
    symmetryY string {mustBeMember(symmetryY, ["PEC", "PMC", "None"])};
    symmetryAxial string {mustBeMember(symmetryAxial, ["TE", "TM", "None"])};
end

%% Get Waveguide Mode Info
modeStructs = nLayer.getCircularModes(...
    O.maxModeIndexM, O.maxModeIndexN, ...
    O.waveguideR, ...
    SymmetryX=symmetryX, SymmetryY=symmetryY, SymmetryAxial=symmetryAxial);

end


