function [] = regenerateWaveguideModeObjects(self)
%Regenerate the "waveguideModes" array.
% The function regeneratates the "waveguideModes" array for an
% "nLayerOpenEnded" object, and will be automatically called in such a way
% to achieve lazy-evaluation.
%
% Author: Matt Dvorsky

arguments
    self nLayerOpenEnded;
end

%% Redefine Mode Structs
if strcmp(class(self), "nLayerOpenEnded")
    self.shouldRegenerateWaveguideModeObjects = false;
    return;
end

if self.shouldRegenerateWaveguideModeObjects
    self.waveguideModes = self.defineWaveguideModes(...
        self.modeSymmetryX, self.modeSymmetryY, self.modeSymmetryAxial);

    self.shouldRecomputeWeights = true;
    self.shouldRegenerateWaveguideModeObjects = false;
end

end

