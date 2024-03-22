function [modeStructs] = defineWaveguideModes(O, symmetryX, symmetryY, symmetryAxial)
%DEFINEWAVEGUIDEMODES Defines waveguide modes for coaxial waveguides.
% Defines the mode spectrums for a coaxial waveguide. Returns a
% modeStruct array as required by the "nLayerOpenEnded" class.
%
% Author: Matt Dvorsky

arguments
    O;
    symmetryX string {mustBeMember(symmetryX, ["PEC", "PMC", "None"])};
    symmetryY string {mustBeMember(symmetryY, ["PEC", "PMC", "None"])};
    symmetryAxial string {mustBeMember(symmetryAxial, ["TE", "TM", "None"])};
end

%% Get Waveguide Mode Info
modeStructs = nLayer.getCoaxialModes(...
    O.maxModeIndexM, O.maxModeIndexN, ...
    O.waveguideRi, O.waveguideRo, ...
    SymmetryX=symmetryX, SymmetryY=symmetryY, SymmetryAxial=symmetryAxial);

end


