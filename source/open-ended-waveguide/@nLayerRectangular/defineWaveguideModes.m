function [waveguideModes] = defineWaveguideModes(O, symmetryX, symmetryY, symmetryAxial)
%Defines waveguide modes for rectangular waveguide.
% Defines the "nLayer.waveguideMode" objects for a rectangular waveguide,
% as required by the "nLayerOpenEnded" class.
%
% Author: Matt Dvorsky

arguments
    O;
    symmetryX string {mustBeMember(symmetryX, ["PEC", "PMC", "None"])};
    symmetryY string {mustBeMember(symmetryY, ["PEC", "PMC", "None"])};
    symmetryAxial string {mustBeMember(symmetryAxial, ["TE", "TM", "None"])};
end

%% Get Waveguide Mode Info
waveguideModes = nLayer.getRectangularModes(...
    O.maxModeIndexM, O.maxModeIndexN, ...
    O.waveguideA, O.waveguideB, ...
    SymmetryX=symmetryX, SymmetryY=symmetryY, ...
    SymmetryAxial=symmetryAxial);

end


