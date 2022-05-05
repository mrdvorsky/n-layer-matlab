function [] = setInitialValues(O, options)
%SETINITIALVALUES Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
    options.ErValue(1, :) {mustBePositive, mustBeFinite, ...
        mustBeGreaterThanOrEqual(options.ErValue, 1)};
    options.ErpValue(1, :) {mustBeNonnegative, mustBeFinite};
    options.ThkValue(1, :) {mustBeNonnegative};
end

%% Check Bounds
if isfield(options, "ThkValue")
    if ~all(isfinite(options.ThkValue(1, 1:end - 1)))
        error("Layer thicknesses must be finite (except for last layer).");
    end
end

%% Assign Initial Guesses
if isfield(options, "ErValue")
    O.erInitialValue =  options.ErValue;
end
if isfield(options, "ErpValue")
    O.erpInitialValue = options.ErpValue;
end
if isfield(options, "ThkValue")
    O.thkInitialValue = options.ThkValue;
end

end

