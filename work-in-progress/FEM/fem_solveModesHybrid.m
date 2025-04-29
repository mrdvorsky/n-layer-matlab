function [valsE, valsH, cutoffs, x, y, tris, isBoundaryPEC, C, T] = fem_solveModesHybrid(x, y, tris, isBoundaryPEC, d, options)
%FEM_SOLVEMODES Summary of this function goes here
% This function ...
%
% Author: Matt Dvorsky

arguments
    x(:, 1) {mustBeReal};
    y(:, 1) {mustBeReal, mustHaveEqualSizes(x, y)};
    tris(:, 3) {mustBeInteger, mustBePositive};
    isBoundaryPEC(:, 1) logical {mustHaveEqualSizes(x, isBoundaryPEC)};
    d(1, 1);

    options.ModeCount(1, 1) {mustBePositive, mustBeInteger} = 1;
    options.ModeType {mustBeMember(options.ModeType, ["TE", "TM", "Hybrid"])} = "TE";

    options.IsBoundaryDir(:, 1) logical {mustHaveEqualSizes(...
        x, options.IsBoundaryDir)} = false(size(x));

    options.FemElementOrder(1, 1) {mustBeInteger, mustBePositive} = 2;

    options.Er(:, 1) = ones(size(tris, 1), 1);
    options.Ur(:, 1) = ones(size(tris, 1), 1);
end

%% Setup Boundary Conditions
if strcmp(options.ModeType, "TM")
    options.IsBoundaryDir = options.IsBoundaryDir | isBoundaryPEC;
end

%% Create Higher Order Triangle Elements
[x, y, tris, isBoundaryPEC] = fem_getGlobalElementCoordinates(...
    options.FemElementOrder, x, y, tris, isBoundaryPEC);
% [x, y, tris] = fem_optimizeNodeConnectivity(x, y, tris);

%% Get Local Matrices
[Tsingle(1, :, :), Qsingle(1, :, :, :), Ysingle(1, :, :, :)] = ...
    fem_generateLocalMatrices(options.FemElementOrder);

%% Create Element-wise Triangle Mesh Descriptors
dy_e = circshift(y(tris(:, 1:3)), -1, 2) - circshift(y(tris(:, 1:3)), -2, 2);
dx_e = circshift(x(tris(:, 1:3)), -2, 2) - circshift(x(tris(:, 1:3)), -1, 2);

% Area of each triangle.
area_e(:, 1) = 0.5 * (dy_e(:, 1).*dx_e(:, 2) - dy_e(:, 2).*dx_e(:, 1));

% Calculate cotangent of each interior angle of each triangle.
cot_e(:, 1, 1, :) = -0.5 * (circshift(dx_e, 1, 2) .* circshift(dx_e, -1, 2) ...
    + circshift(dy_e, 1, 2) .* circshift(dy_e, -1, 2)) ...
    ./ area_e;

%% Construct Matrices C and T Element-wise
% Calculate element matrices C and T. See Sadiku p393, p413, p416. Note
% that on p416, there is a mistake in (6.113.5). The denominator of 2A
% should not be present.
C_e = sum(Qsingle .* cot_e, 4);
Y_e = sum(Ysingle .* cot_e, 4);
T_e = area_e .* Tsingle;

%% Construct Global A and B from Element-wise A and B
ind_i = tris + 0*reshape(tris, size(tris, 1), 1, []);
ind_j = 0*tris + reshape(tris, size(tris, 1), 1, []);

kt2 = options.Er .* options.Ur - d.^2;

C_E = sparse(ind_i(:), ind_j(:), reshape(options.Er ./ kt2 .* C_e, [], 1));
T_E = sparse(ind_i(:), ind_j(:), reshape(options.Er .* T_e, [], 1));

C_H = sparse(ind_i(:), ind_j(:), reshape(options.Ur ./ kt2 .* C_e, [], 1));
T_H = sparse(ind_i(:), ind_j(:), reshape(options.Ur .* T_e, [], 1));

Y = sparse(ind_i(:), ind_j(:), reshape(d ./ kt2 .* Y_e, [], 1));
Y_E = Y;
Y_H = Y;
Y_E(isBoundaryPEC, :) = [];
Y_H(:, isBoundaryPEC) = [];

C_E(isBoundaryPEC, :) = [];
T_E(isBoundaryPEC, :) = [];
C_E(:, isBoundaryPEC) = [];
T_E(:, isBoundaryPEC) = [];

s = [size(C_E, 1), size(C_H, 2)];
C = [C_E, Y_E; Y_H, C_H];
T = [T_E, sparse(s(1), s(2)); sparse(s(2), s(1)), T_H];

[V, D] = eigs(C, T, options.ModeCount, ...
    "smallestabs", IsSymmetricDefinite=true, Display=true);

valsE = zeros(numel(x), size(V, 2));
valsE(~isBoundaryPEC, :) = V(1:s(1), :);
valsH = V(s(1) + 1:end, :);

cutoffs = sqrt(diag(D));

end

