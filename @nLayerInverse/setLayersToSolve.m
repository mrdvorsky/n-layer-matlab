function [] = setLayersToSolve(O, options)
%SETLAYERSTOSOLVE Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
    options.ErLayers(:, 1) {mustBeInteger, mustBePositive} = [];
    options.ErpLayers(:, 1) {mustBeInteger, mustBePositive} = [];
    options.ThkLayers(:, 1) {mustBeInteger, mustBePositive} = [];
end

%% Check Index Range
if any(options.ErLayers > O.layerCount)
    error(strcat("All elements of 'ErLayers' must be not greater ", ...
        "than the number of layers (%d)."), O.layerCount);
end

if any(options.ErpLayers > O.layerCount)
    error(strcat("All elements of 'ErpLayers' must be not greater ", ...
        "than the number of layers (%d)."), O.layerCount);
end

if any(options.ThkLayers > O.layerCount)
    error(strcat("All elements of 'ThkLayers' must be not greater ", ...
        "than the number of layers (%d)."), O.layerCount);
end

%% Assign Layer Indices
O.erLayersToSolve =  options.ErLayers;
O.erpLayersToSolve = options.ErpLayers;
O.thkLayersToSolve = options.ThkLayers;

end

