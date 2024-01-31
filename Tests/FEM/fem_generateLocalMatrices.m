function [T, Q, Y] = fem_generateLocalMatrices(orderNum)
%FEM_GENERATELOCALMATRICES Generate local matrices T and Q for 2D FEM.
% This function ...
%
% Author: Matt Dvorsky

arguments
    orderNum(1, 1) {mustBeInteger, mustBePositive};
end

%% Get Shape Functions
shapeFun = fem_generateShapeFunctions(orderNum);

%% Get Weights and Nodes for Triangle
[nodes, weights] = gaussLegendre(orderNum + 1, 0, 1);

nodesX1 = reshape(nodes + 0*nodes.', [], 1);
nodesX2 = reshape(nodes.' .* (1 - nodes), [], 1);
weightsX12 = reshape(weights .* weights.' .* (1 - nodes), [], 1);

%% Evaluate T
% Evaluate
TSamp3 = fem_evaluateShapeFunctions(shapeFun, nodesX1, nodesX2);
TSamp4 = reshape(TSamp3, size(TSamp3, 1), 1, 1, []);

T = shiftdim(innerProduct(...
    2 * weightsX12 .* TSamp3, ...
    TSamp4, ...
    1), 2);

% Force symmetry to eliminate numerical error.
T = 0.5 * (T + T.');

%% Calculate Derivatives of ShapeFun
shapeFunDerAll = [shapeFun(:, 2:end, :) .* (1:size(shapeFun, 2) - 1), ...
    0*shapeFun(:, 1, :)];

shapeFunDer1 = [shapeFunDerAll(1, :, :); shapeFun(2, :, :); shapeFun(3, :, :)];
shapeFunDer2 = [shapeFun(1, :, :); shapeFunDerAll(2, :, :); shapeFun(3, :, :)];
shapeFunDer3 = [shapeFun(1, :, :); shapeFun(2, :, :); shapeFunDerAll(3, :, :)];

%% Calculate Q
QSampDer1_3 = fem_evaluateShapeFunctions(shapeFunDer1, nodesX1, nodesX2);
QSampDer2_3 = fem_evaluateShapeFunctions(shapeFunDer2, nodesX1, nodesX2);
QSampDer3_3 = fem_evaluateShapeFunctions(shapeFunDer3, nodesX1, nodesX2);
QSampDer1_4 = reshape(QSampDer1_3, size(QSampDer1_3, 1), 1, 1, []);
QSampDer2_4 = reshape(QSampDer2_3, size(QSampDer2_3, 1), 1, 1, []);
QSampDer3_4 = reshape(QSampDer3_3, size(QSampDer3_3, 1), 1, 1, []);

Q1 = shiftdim(innerProduct(...
    weightsX12 .* (QSampDer2_3 - QSampDer3_3), ...
    (QSampDer2_4 - QSampDer3_4), ...
    1), 2);
Q2 = shiftdim(innerProduct(...
    weightsX12 .* (QSampDer3_3 - QSampDer1_3), ...
    (QSampDer3_4 - QSampDer1_4), ...
    1), 2);
Q3 = shiftdim(innerProduct(...
    weightsX12 .* (QSampDer1_3 - QSampDer2_3), ...
    (QSampDer1_4 - QSampDer2_4), ...
    1), 2);

Q = cat(3, Q1, Q2, Q3);

% Force symmetry to eliminate numerical error.
Q = 0.5 * (Q + pagetranspose(Q));

%% Calculate Y
Y1 = shiftdim(innerProduct(...
    weightsX12 .* (QSampDer2_3 + QSampDer3_3), ...
    (QSampDer2_4 + QSampDer3_4), ...
    1), 2);
Y2 = shiftdim(innerProduct(...
    weightsX12 .* (QSampDer3_3 + QSampDer1_3), ...
    (QSampDer3_4 + QSampDer1_4), ...
    1), 2);
Y3 = shiftdim(innerProduct(...
    weightsX12 .* (QSampDer1_3 + QSampDer2_3), ...
    (QSampDer1_4 + QSampDer2_4), ...
    1), 2);

Y = cat(3, Y1, Y2, Y3);

% Force symmetry to eliminate numerical error.
Y = 0.5 * (Y + pagetranspose(Y));

end

