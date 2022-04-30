function [] = setLayerCount(O, layerCount)
%SETLAYERCOUNT Sets the layer count for the multilayer structure.
%   Detailed explanation goes here

arguments
    O;
    layerCount(1, 1) {mustBeInteger, mustBePositive};
end

%% Check Default Parameters


%% Set Ranges and Guesses
O.layerCount = layerCount;

O.erRange = O.default_erRange + 0*(1:layerCount);
O.erpRange = O.default_erpRange + 0*(1:layerCount);
O.thkRange = O.default_thkRange + 0*(1:layerCount);

O.erGuess = O.erRange(1, :);
O.erpGuess = O.erpRange(1, :);
O.thkGuess = O.thkRange(1, :);

end

