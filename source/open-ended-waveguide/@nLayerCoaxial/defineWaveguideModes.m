function [modes] = defineWaveguideModes(self, symmetryX, symmetryY, symmetryAxial)
%Defines waveguide modes for a coaxial waveguide.
% Defines the "nLayer.waveguideMode" objects for a coaxial waveguide,
% as required by the "nLayerOpenEnded" class.
%
% Author: Matt Dvorsky

arguments
    self nLayerCoaxial;

    symmetryX string {mustBeMember(symmetryX, ["PEC", "PMC", "None"])};
    symmetryY string {mustBeMember(symmetryY, ["PEC", "PMC", "None"])};
    symmetryAxial string {mustBeMember(symmetryAxial, ["TE", "TM", "None"])};
end

%% Get Waveguide Mode Info
modes = waveguideMode.getAllCoaxialModes(...
    0:self.modeIndexM, 0:self.maxModeIndexN, ...
    self.waveguideRi, self.waveguideRo, ...
    SymmetryX=symmetryX, SymmetryY=symmetryY, ...
    SymmetryAxial=symmetryAxial);

end


