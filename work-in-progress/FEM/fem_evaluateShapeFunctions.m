function [y] = fem_evaluateShapeFunctions(shapeFun, xi1, xi2, xi3)
%FEM_EVALUATESHAPEFUNCTIONS Evaluate 2D FEM multinomial shape functions.
% This function ...
%
% Author: Matt Dvorsky

arguments
    shapeFun {mustBeValidShapeFun};
    xi1 {mustBeInRange(xi1, 0, 1), ...
        mustHaveCompatibleSizes(shapeFun, xi1, ExcludeDimensions=1:2)};
    xi2 {mustBeInRange(xi2, 0, 1), mustHaveCompatibleSizes(xi1, xi2), ...
        mustHaveCompatibleSizes(shapeFun, xi1, xi2, ExcludeDimensions=1:2)};
    xi3 {mustBeInRange(xi3, 0, 1), mustHaveCompatibleSizes(xi1, xi2, xi3), ...
        mustHaveCompatibleSizes(shapeFun, xi1, xi2, xi3, ExcludeDimensions=1:2)} ...
        = 1 - xi1 - xi2;
end

%% Check Inputs and Shift Dimensions
xi1 = shiftdim(xi1, -2);
xi2 = shiftdim(xi2, -2);
xi3 = shiftdim(xi3, -2);

shapeFunSizeOld = size(shapeFun);
shapeFunSize = [shapeFunSizeOld(1:2), 1, 1, shapeFunSizeOld(3:end)];
shapeFun = reshape(shapeFun, shapeFunSize);

%% Evaluate Multinomial
powInd(1, :) = (0:size(shapeFun, 2) - 1);
indAll = repmat({':'}, 1, ndims(shapeFun) - 1);

y = shiftdim(...
       innerProduct(shapeFun(1, indAll{:}), xi1.^powInd, 2) ...
    .* innerProduct(shapeFun(2, indAll{:}), xi2.^powInd, 2) ...
    .* innerProduct(shapeFun(3, indAll{:}), xi3.^powInd, 2), ...
    2);

end




%% Custom Argument Validation Function
function mustBeValidShapeFun(shapeFun)
    if size(shapeFun, 1) ~= 3
        throwAsCaller(MException("MATLAB:mustBeValidShapeFun", ...
            "First argument must be of size 3-by-n-by-..."));
    end
end

