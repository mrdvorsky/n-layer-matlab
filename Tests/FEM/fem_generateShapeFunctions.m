function [shapeMultinomials] = fem_generateShapeFunctions(orderNum)
%FEM_CREATESHAPEFUNCTIONS Create shape functions for 2D FEM solver.
% This function ...
%
% Author: Matt Dvorsky

arguments
    orderNum(1, 1) {mustBeInteger, mustBePositive};
end

%% Create IJK Index List in Sorted Order (Sadiku p413)
% Create in the order defined by Sadiku p413.
[indI_tmp, indJ_tmp] = ndgrid(0:orderNum, 0:orderNum);
indSIJ = [indI_tmp(:) + indJ_tmp(:), indI_tmp(:), indJ_tmp(:)];
indSIJ = sortrows(indSIJ);
indIJK = [orderNum - indSIJ(:, 1), indSIJ(:, 3), indSIJ(:, 2)];
indIJK = indIJK(indIJK(:, 1) >= 0, :);

% Sort so that the order is vertices, edges, interior.
verticesInd = [1, 0.5*(orderNum + 1)*(orderNum + 2) - [orderNum, 0]];
edgesInd = [...
    cumsum(1:orderNum - 1) + 1, ...
    (verticesInd(2) + 1):(verticesInd(3) - 1), ...
    flip(cumsum(2:orderNum) + 1)];
interiorInd = find(all(indIJK(:, :) ~= 0, 2)).';

indIJK = indIJK([verticesInd, edgesInd, interiorInd], :);

%% Create Shape Functions Helpers (Sadiku p412)
pr = cell(orderNum + 1, 1);
pr{1} = 1 * ((0:orderNum) == 0);
for r = 1:orderNum
    pr{r + 1} = resize(conv(pr{r}, [1 - r, orderNum] ./ r), orderNum + 1);
end

%% Create Shape Functions
% Shape functions are multinomials in xi1, xi2, xi3 (e.g., Sadiku p414).
%  shapeMultinomials(j, :, n) is the polynomial in xi_j for the nth shape
%  function. Full shape function is a product of the 3 polynomials.
shapeMultinomials = zeros(3, orderNum + 1, size(indIJK, 1));
for ii = 1:size(indIJK, 1)
    shapeMultinomials(1, :, ii) = pr{indIJK(ii, 1) + 1};
    shapeMultinomials(2, :, ii) = pr{indIJK(ii, 2) + 1};
    shapeMultinomials(3, :, ii) = pr{indIJK(ii, 3) + 1};
end


end





