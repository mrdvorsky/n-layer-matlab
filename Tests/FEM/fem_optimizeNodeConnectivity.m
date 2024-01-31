function [xOut, yOut, trisOut] = fem_optimizeNodeConnectivity(x, y, tris)
%FEM_OPTIMIZENODECONNECTIVITY Permute nodes to get optimal sparse matrix.
% This function ...
%
% Author: Matt Dvorsky

arguments
    x(:, 1) {mustBeReal};
    y(:, 1) {mustBeReal, mustHaveEqualSizes(x, y)};
    tris(:, :) {mustBeInteger, mustBePositive};
end

%% Build Sparse Matrix and Find Optimal Permutation Vector
ind_i = tris + 0*reshape(tris, size(tris, 1), 1, []);
ind_j = 0*tris + reshape(tris, size(tris, 1), 1, []);
permInd = symrcm(sparse(ind_i(:), ind_j(:), ones(size(ind_i(:)))));

%% Apply Permutation to Inputs
xOut = x(permInd);
yOut = y(permInd);

permIndInv(permInd) = 1:numel(permInd);
trisOut = reshape(permIndInv(tris), size(tris));

end

