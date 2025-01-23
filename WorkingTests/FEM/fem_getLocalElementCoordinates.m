function [xi1, xi2, xi3] = fem_getLocalElementCoordinates(orderNum)
%FEM_GETLOCALELEMENTCOORDINATES Get local coordinates of 
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

%% Set Outputs
xi123 = indIJK ./ orderNum;

if nargout <= 1
    xi1 = xi123;
else
    xi1 = xi123(:, 1);
    xi2 = xi123(:, 2);
    xi3 = xi123(:, 3);
end

end

