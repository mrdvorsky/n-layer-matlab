function [] = regenerateModeStructs(O)
%REGENERATEMODESTRUCTS Regenerate modeStructs array.
% The function regeneratates the "modeStructs" array for an
% "nLayerOpenEnded" object, and will be automatically called whenever a
% parameter changes that (...).
%
% Inputs:
%   parameterName - A string containing the name of the parameter that
%       changed. This can be used to customize behavior using a switch
%       statement for optimization for other reasons.
%
% Author: Matt Dvorsky

arguments
    O;
end

%% Redefine Mode Structs
O.modeStructs = O.defineWaveguideModes(...
    O.modeSymmetryX, O.modeSymmetryY, O.modeSymmetryAxial);

end

