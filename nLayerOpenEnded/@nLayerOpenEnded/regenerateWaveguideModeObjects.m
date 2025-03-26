function [] = regenerateWaveguideModeObjects(O)
%REGENERATEWAVEGUIDEMODEOBJECTS Regenerate the "waveguideModes" array.
% The function regeneratates the "waveguideModes" array for an
% "nLayerOpenEnded" object, and will be automatically called in such a way
% to achieve lazy-evaluation.
%
% Author: Matt Dvorsky

arguments
    O;
end

%% Redefine Mode Structs
if strcmp(class(O), "nLayerOpenEnded")
    O.shouldRegenerateWaveguideModeObjects = false;
    return;
end

if O.shouldRegenerateWaveguideModeObjects
    O.waveguideModes = O.defineWaveguideModes(...
        O.modeSymmetryX, O.modeSymmetryY, O.modeSymmetryAxial);

    O.shouldRecomputeWeights = true;
    O.shouldRegenerateWaveguideModeObjects = false;
end

end

