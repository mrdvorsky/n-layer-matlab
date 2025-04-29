function [waveguideModes] = defineWaveguideModes(O, symmetryX, symmetryY, symmetryAxial)
%Defines waveguide modes for a circular waveguide.
% Defines the "nLayer.waveguideMode" objects for a circular waveguide,
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
waveguideModes = nLayer.getCircularModes(...
    O.modeIndexM, O.maxModeIndexN, ...
    O.waveguideR, ...
    SymmetryX=symmetryX, SymmetryY=symmetryY, ...
    SymmetryAxial=symmetryAxial);

end


