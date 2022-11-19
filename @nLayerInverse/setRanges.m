function [] = setRanges(O, options)
%SETRANGES Sets the min and max ranges for each variable and layer.
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
    O;
    options.ErpMin (:, 1) {mustBeValidRangeMin(options.ErpMin,  O)} = O.rangeMin_erp;
    options.ErppMin(:, 1) {mustBeValidRangeMin(options.ErppMin, O)} = O.rangeMin_erpp;
    options.UrpMin (:, 1) {mustBeValidRangeMin(options.UrpMin,  O)} = O.rangeMin_urp;
    options.UrppMin(:, 1) {mustBeValidRangeMin(options.UrppMin, O)} = O.rangeMin_urpp;
    options.ThkMin (:, 1) {mustBeValidRangeMin(options.ThkMin,  O)} = O.rangeMin_thk;
    
    options.ErpMax (:, 1) {mustBeValidRangeMax(options.ErpMax,  O)} = O.rangeMax_erp;
    options.ErppMax(:, 1) {mustBeValidRangeMax(options.ErppMax, O)} = O.rangeMax_erpp;
    options.UrpMax (:, 1) {mustBeValidRangeMax(options.UrpMax,  O)} = O.rangeMax_urp;
    options.UrppMax(:, 1) {mustBeValidRangeMax(options.UrppMax, O)} = O.rangeMax_urpp;
    options.ThkMax (:, 1) {mustBeValidRangeMax(options.ThkMax,  O)} = O.rangeMax_thk;
end

%% Assign Ranges
O.rangeMin_erp  = options.ErpMin;
O.rangeMin_erpp = options.ErppMin;
O.rangeMin_urp  = options.UrpMin;
O.rangeMin_urpp = options.UrppMin;
O.rangeMin_thk  = options.ThkMin;

O.rangeMax_erp  = options.ErpMax;
O.rangeMax_erpp = options.ErppMax;
O.rangeMax_urp  = options.UrpMax;
O.rangeMax_urpp = options.UrppMax;
O.rangeMax_thk  = options.ThkMax;

end

function mustBeValidRangeMin(ranges, O)
    if numel(ranges) ~= O.layerCount
        throwAsCaller(MException("nLayerInverse:mustBeValidRangeMin", ...
            "The parameters 'rangeMin_{er, ur, thk}' must be vectors " + ...
            "with layerCount (%d) elements.", O.layerCount));
    end
end

function mustBeValidRangeMax(ranges, O)
    if numel(ranges) ~= O.layerCount
        throwAsCaller(MException("nLayerInverse:mustBeValidRangeMax", ...
            "The parameters 'rangeMax_{er, ur, thk}' must be vectors " + ...
            "with layerCount (%d) elements.", O.layerCount));
    end
end

