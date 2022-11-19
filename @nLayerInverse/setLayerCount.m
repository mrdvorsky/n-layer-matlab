function [] = setLayerCount(O, layerCount)
%SETLAYERCOUNT Changes the layer count for the multilayer structure.
% Changes the layerCount property. Additionally, modifies the initial
% value and range parameters to be the correct size. If new layers are
% being added, will duplicate the current last layer values. Otherwise,
% will not change layer values.
%
% Example Usage:
%   NLsolver.setLayerCount(numLayers);
%
% Author: Matt Dvorsky

arguments
    O;
    layerCount(1, 1) {mustBeInteger, mustBePositive};
end

%% Set 'layerCount' Property
O.layerCount = layerCount;

%% Change Initial Value and Range Parameters to be the Correct Size
O.rangeMin_erp  = updateValues(O.rangeMin_erp,  layerCount);
O.rangeMin_erpp = updateValues(O.rangeMin_erpp, layerCount);
O.rangeMin_urp  = updateValues(O.rangeMin_urp,  layerCount);
O.rangeMin_urpp = updateValues(O.rangeMin_urpp, layerCount);
O.rangeMin_thk  = updateValues(O.rangeMin_thk,  layerCount);

O.rangeMax_erp  = updateValues(O.rangeMax_erp,  layerCount);
O.rangeMax_erpp = updateValues(O.rangeMax_erpp, layerCount);
O.rangeMax_urp  = updateValues(O.rangeMax_urp,  layerCount);
O.rangeMax_urpp = updateValues(O.rangeMax_urpp, layerCount);
O.rangeMax_thk  = updateValues(O.rangeMax_thk,  layerCount);

O.initialValue_er  = updateValues(O.initialValue_er,  layerCount);
O.initialValue_ur  = updateValues(O.initialValue_ur,  layerCount);
O.initialValue_thk = updateValues(O.initialValue_thk, layerCount);

end

function [valuesOut] = updateValues(valuesIn, layerCount)
    if numel(valuesIn) >= layerCount
        valuesOut = valuesIn(1, 1:layerCount);
    else
        valuesOut = [valuesIn, ...
            repmat(valuesIn(1, end), 1, layerCount - numel(valuesIn))];
    end
end


