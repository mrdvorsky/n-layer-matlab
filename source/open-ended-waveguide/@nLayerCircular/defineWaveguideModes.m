function [modes] = defineWaveguideModes(self, symmetryX, symmetryY, symmetryAxial)
%Defines waveguide modes for a circular waveguide.
% Defines the "waveguideMode" objects for a circular waveguide,
% as required by the "nLayerOpenEnded" class.
%
% Author: Matt Dvorsky

arguments
    self nLayerCircular;

    symmetryX string {mustBeMember(symmetryX, ["PEC", "PMC", "None"])};
    symmetryY string {mustBeMember(symmetryY, ["PEC", "PMC", "None"])};
    symmetryAxial string {mustBeMember(symmetryAxial, ["TE", "TM", "None"])};
end

%% Get Waveguide Mode Info
modes = waveguideMode.getAllCircularModes(...
    0:self.modeIndexM, 1:self.maxModeIndexN, ...
    self.waveguideR, ...
    SymmetryX=symmetryX, SymmetryY=symmetryY, ...
    SymmetryAxial=symmetryAxial);

end


