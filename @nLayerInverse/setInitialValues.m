function [] = setInitialValues(O, options)
%SETINITIALVALUES Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
    options.ErValue(1, :) {mustBePositive, mustBeFinite, ...
        mustBeGreaterThanOrEqual(options.ErValue, 1)};
    options.ErpValue(1, :) {mustBePositive, mustBeFinite};
    options.ThkValue(1, :) {mustBePositive};
end

%% Check Bounds
if any(~isfinite(options.ThkValue(1, 1:end - 1)))
    error("Layer thicknesses must be finite (except for last layer).");
end

%% Assign Initial Guesses
O.erInitialValue =  options.ErValue;
O.erpInitialValue = options.ErpValue;
O.thkInitialValue = options.ThkValue;

end

