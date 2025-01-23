function [xOut, yOut, trisOut, varargout] = fem_getGlobalElementCoordinates(orderNum, x, y, tris, isBoundary)
%FEM_GETGLOBALELEMENTCOORDINATES Get coordinates of shape function nodes.
% This function ...
%
% Author: Matt Dvorsky

arguments
    orderNum(1, 1) {mustBeInteger, mustBePositive};
    x(:, 1) {mustBeReal};
    y(:, 1) {mustBeReal, mustHaveEqualSizes(x, y)};
    tris(:, 3) {mustBeInteger, mustBePositive};
end
arguments (Repeating)
    isBoundary(:, 1) logical;
end

%% Check Inputs
if orderNum == 1
    xOut = x;
    yOut = y;
    trisOut = tris;
    varargout = isBoundary;
    return;
end

%% Compute Edge Connectivity Matrix
[edges, indTmp] = sort([tris(:), reshape(circshift(tris, -1, 2), [], 1)], 2);
isEdgeFlipped = reshape(indTmp(:, 1) ~= 1, [], 1);

[edgesSorted, edgeInd] = sortrows(edges);

[~, ia, ~] = unique(edgesSorted, "first", "rows");
dupInds = setdiff(1:size(edgesSorted, 1), ia);

edgeNodes = zeros(size(edgesSorted, 1), orderNum - 1);
edgeNodes(ia, :) = max(edgesSorted, [], "all") ...
    + reshape(1:numel(ia)*(orderNum - 1), orderNum - 1, []).';
edgeNodes(dupInds, :) = edgeNodes(dupInds - 1, :);
dupIndsFlip = xor(isEdgeFlipped(edgeInd(dupInds)), isEdgeFlipped(edgeInd(dupInds - 1)));
edgeNodes(dupInds(dupIndsFlip), :) = flip(edgeNodes(dupInds(dupIndsFlip), :), 2);

for ii = 1:numel(isBoundary)
    isEdgeBoundary = all(isBoundary{ii}(edgesSorted), 2);
    isEdgeBoundary([dupInds - 1, dupInds]) = false;
    boundaryEdgeInds{ii} = unique(edgeNodes(isEdgeBoundary, :));
end

edgesNodesNew(edgeInd, :) = edgeNodes;
edgesNodesNew = permute(reshape(edgesNodesNew, [], 3, orderNum - 1), [1, 3, 2]);

interiorNodes = zeros(0.5 * (orderNum - 1)*(orderNum - 2), size(tris, 1));
interiorNodes(:) = 1:numel(interiorNodes);
interiorNodes = max(edgesNodesNew, [], "all") + interiorNodes.';

trisOut = [tris, edgesNodesNew(:, :), interiorNodes];

%% Calculate Coordinates of Higher Order Nodes
x123 = reshape(x(tris), [], 1, 3);
y123 = reshape(y(tris), [], 1, 3);

xi123(1, :, :) = fem_getLocalElementCoordinates(orderNum);

xAll = innerProduct(xi123, x123, 3);
yAll = innerProduct(xi123, y123, 3);

%% Calculate Output Node Locations
xOut(trisOut(:), 1) = xAll(:);
yOut(trisOut(:), 1) = yAll(:);

%% Compute Boundary Specifiers
varargout = cell(numel(isBoundary), 1);
for ii = 1:numel(isBoundary)
    varargout{ii} = resize(isBoundary{ii}, numel(xOut));
    varargout{ii}(boundaryEdgeInds{ii}) = true;
end

end

