function [] = regenerateModeStructs(O)
%REGENERATEMODESTRUCTS Regenerate modeStructs array.
% The function regeneratates the "modeStructs" array for an
% "nLayerOpenEnded" object, and will be automatically called whenever a
% parameter changes that (...).
%
% Author: Matt Dvorsky

arguments
    O;
end

%% Redefine Mode Structs
wasEmpty = true;
if ~isempty(O.modeStructs)
    er = O.waveguideEr;
    ur = O.waveguideUr;
    wasEmpty = false;
end

O.modeStructs = O.defineWaveguideModes(...
    O.modeSymmetryX, O.modeSymmetryY, O.modeSymmetryAxial);

if ~wasEmpty
    O.waveguideEr = er;
    O.waveguideUr = ur;
end

end

