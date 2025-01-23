function [vals, cutoffs, x, y, tris, isBoundaryPEC, C, T] = fem_solveModes(x, y, tris, isBoundaryPEC, options)
%FEM_SOLVEMODES Summary of this function goes here
% This function ...
%
% Author: Matt Dvorsky

arguments
    x(:, 1) {mustBeReal};
    y(:, 1) {mustBeReal, mustHaveEqualSizes(x, y)};
    tris(:, 3) {mustBeInteger, mustBePositive};
    isBoundaryPEC(:, 1) logical {mustHaveEqualSizes(x, isBoundaryPEC)};

    options.ModeCount(1, 1) {mustBePositive, mustBeInteger} = 1;
    options.ModeType {mustBeMember(options.ModeType, ["TE", "TM", "Hybrid"])} = "TE";

    options.IsBoundaryDir(:, 1) logical {mustHaveEqualSizes(...
        x, options.IsBoundaryDir)} = false(size(x));

    options.FemElementOrder(1, 1) {mustBeInteger, mustBePositive} = 2;

    options.Er(:, 1) {mustHaveEqualSizes(options.Er, x)} = ones(size(x));
    options.Ur(:, 1) {mustHaveEqualSizes(options.Ur, x)} = ones(size(x));
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
[Tsingle(1, :, :), Qsingle(1, :, :, :)] = ...
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
T_e = area_e .* Tsingle;

if strcmp(options.ModeType, "TM")
    % Construct Global A and B from Element-wise A and B
    ind_i = tris + 0*reshape(tris, size(tris, 1), 1, []);
    ind_j = 0*tris + reshape(tris, size(tris, 1), 1, []);

    C = sparse(ind_i(:), ind_j(:), C_e(:));
    T = sparse(ind_i(:), ind_j(:), T_e(:));

    % Solve Global System
    [V, D] = eigs(C, T, options.ModeCount + 1, ...
        "smallestabs", IsSymmetricDefinite=true);
    vals = V(:, 2:end);
    cutoffs = sqrt(diag(D(2:end, 2:end)));
elseif strcmp(options.ModeType, "TE")
    % Construct Global A and B from Element-wise A and B
    ind_i = tris + 0*reshape(tris, size(tris, 1), 1, []);
    ind_j = 0*tris + reshape(tris, size(tris, 1), 1, []);

    C = sparse(ind_i(:), ind_j(:), C_e(:));
    T = sparse(ind_i(:), ind_j(:), T_e(:));

    C(isBoundaryPEC, :) = [];
    T(isBoundaryPEC, :) = [];

    C(:, isBoundaryPEC) = [];
    T(:, isBoundaryPEC) = [];

    % Solve Global System
    [V, D] = eigs(C, T, options.ModeCount + 1, ...
        "smallestabs", IsSymmetricDefinite=true);

    vals = zeros(numel(isBoundaryPEC), size(V, 2));
    vals(~isBoundaryPEC, :) = V;
    cutoffs = sqrt(diag(D));
else
    % Construct Global A and B from Element-wise A and B
    ind_i = tris + 0*reshape(tris, size(tris, 1), 1, []);
    ind_j = 0*tris + reshape(tris, size(tris, 1), 1, []);

    C_TE = sparse(ind_i(:), ind_j(:), reshape(C_e, [], 1));
    T_TE = sparse(ind_i(:), ind_j(:), reshape(T_e, [], 1));

    C_TM = sparse(ind_i(:), ind_j(:), reshape(C_e, [], 1));
    T_TM = sparse(ind_i(:), ind_j(:), reshape(T_e, [], 1));

    C_TE(isBoundaryPEC, :) = [];
    T_TE(isBoundaryPEC, :) = [];
    C_TE(:, isBoundaryPEC) = [];
    T_TE(:, isBoundaryPEC) = [];

    s = [size(C_TE, 1), size(C_TM, 2)];
    C = [C_TE, sparse(s(1), s(2)); sparse(s(2), s(1)), C_TM];
    T = [T_TE, sparse(s(1), s(2)); sparse(s(2), s(1)), T_TM];

    [V, D] = eigs(C, T, options.ModeCount + 1, ...
        "smallestabs", IsSymmetricDefinite=true);
end

end

