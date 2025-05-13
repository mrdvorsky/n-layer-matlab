function [waveguideModes] = defineWaveguideModes(self, symmetryX, symmetryY, symmetryAxial)
%Defines waveguide modes for rectangular waveguide.
% Defines the "nLayer.waveguideMode" objects for a rectangular waveguide,
% as required by the "nLayerOpenEnded" class.
%
% Author: Matt Dvorsky

arguments
    self nLayerRectangular;

    symmetryX string {mustBeMember(symmetryX, ["PEC", "PMC", "None"])};
    symmetryY string {mustBeMember(symmetryY, ["PEC", "PMC", "None"])};
    symmetryAxial string {mustBeMember(symmetryAxial, ["TE", "TM", "None"])};
end

%% Get Waveguide Mode Info
waveguideModes = waveguideMode.getAllRectangularModes(...
    0:self.maxModeIndexM, 0:self.maxModeIndexN, ...
    self.waveguideA, self.waveguideB, ...
    SymmetryX=symmetryX, SymmetryY=symmetryY, ...
    SymmetryAxial=symmetryAxial);

end


