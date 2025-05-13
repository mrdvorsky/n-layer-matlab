function [] = setRanges(self, options)
%Sets the min and max ranges for each variable and layer.
% This function is simply a setter for the rangeMin_{...} and
% rangeMax_{...} parameters that does some bounds and size checking. Any
% unspecified parameter ranges will be left unchanged.
%
% Example Usage (for 2 layers):
%   NLsolver.setRanges(ErpMin=[1, 5], ErpMax=[inf, 10]);
%   NLsolver.setRanges(ThkMax=[10, inf], ErppMax=[0.01, 10]);
%
% Author: Matt Dvorsky

arguments
    self nLayerInverse;

    options.ErpMin (1, :) = self.rangeMin_erp;
    options.ErppMin(1, :) = self.rangeMin_erpp;
    options.UrpMin (1, :) = self.rangeMin_urp;
    options.UrppMin(1, :) = self.rangeMin_urpp;
    options.ThkMin (1, :) = self.rangeMin_thk;

    options.ErpMax (1, :) = self.rangeMax_erp;
    options.ErppMax(1, :) = self.rangeMax_erpp;
    options.UrpMax (1, :) = self.rangeMax_urp;
    options.UrppMax(1, :) = self.rangeMax_urpp;
    options.ThkMax (1, :) = self.rangeMax_thk;
end

%% Assign Ranges
self.rangeMin_erp  = options.ErpMin;
self.rangeMin_erpp = options.ErppMin;
self.rangeMin_urp  = options.UrpMin;
self.rangeMin_urpp = options.UrppMin;
self.rangeMin_thk  = options.ThkMin;

self.rangeMax_erp  = options.ErpMax;
self.rangeMax_erpp = options.ErppMax;
self.rangeMax_urp  = options.UrpMax;
self.rangeMax_urpp = options.UrppMax;
self.rangeMax_thk  = options.ThkMax;

end

