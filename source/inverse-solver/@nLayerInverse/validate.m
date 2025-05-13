function [] = validate(self)
%Checks that the nLayerInverse object is valid.
% This function checks that the sizes of the initialValue_{...},
% rangeMin_{...}, and rangeMax_{...} are vectors with a length of
% layerCount. Also checks that the layersToSolve_{...} parameters have
% unique integers no larger than layerCount.
%
% Since the 'setInitialValues', 'setLayersToSolve', and 'setRanges'
% functions already implement this checking, this function simply calls
% each with no arguments.
%
% Example Usage:
%   NLsolver.validate();
%
%
% Author: Matt Dvorsky

arguments
    self nLayerInverse;
end

%% Check Validity
% Use the setter functions, since they already check these things.
self.setInitialValues();
self.setLayersToSolve();
self.setRanges();

end

