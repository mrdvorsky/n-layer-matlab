function [] = addThicknessConstraint(O, layerIndices, constraints)
%ADDTHICKNESSCONSTRAINT Add thickness constraints to structure.
% This functions adds thickness constraints to an "nLayerInverse" solver.
% Typically, this function is used to constrain the total structure
% thickness, but it can be used to constrain any number of layers with
% either equality or comparison constraints.
%
% Example Usage (for plug with unknown thickness):
%   NLsolver.setInitialValues(Er=[1, 4 - 0.1j, 1], Thk=[5, 5, 5]);
%   NLsolver.setLayersToSolve(Er=[2], Thk=[1, 2, 3]);
%
%   % Next line enforces that total structure thickness won't change.
%   NLsolver.addThicknessConstraint("all", IsFixed=true);
%   
%   % Other examples.
%   NLsolver.addThicknessConstraint([1, 2, 3], LessThanOrEqualTo=15);
%   NLsolver.addThicknessConstraint([2, 3],    GreaterThanOrEqualTo=5);
%
%
% Named Arguments:
%   SumOfLayers ("all") - Layer indices of thicknesses to sum.
%   IsFixed (false) - If true, the total thickness of the specified layers
%       will not change during optimization.
%   LessThanOrEqualTo - Scalar thickness that the sum should be less than
%       or equal to. Only specify one of these scalar constraints.
%   GreaterThanOrEqualTo - Scalar thickness that the sum should be greater
%       than or equal to. Only specify one of these scalar constraints.
%
% Author: Matt Dvorsky

arguments
    O;
    layerIndices(1, :) {mustBeValidDimension};
    constraints.IsFixed(1, 1) logical = false;
    constraints.LessThanOrEqualTo(1, 1) {mustBePositive};
    constraints.GreaterThanOrEqualTo(1, 1) {mustBePositive};
end

%% Check Inputs
if strcmp(layerIndices, "all")
    layerIndices = 1:O.layerCount;
end

if max(layerIndices) > O.layerCount
    error("Input 'layerIndices' must be a vector of valid layer " + ...
        "indices or the string 'all'.");
end

%% Set Constraint Matrices
thk_A = zeros(1, O.layerCount);
thk_A(layerIndices) = 1;

if constraints.IsFixed
    O.constraints_thk_Aeq = [O.constraints_thk_Aeq; thk_A];
    return;
end

if isfield(constraints, "LessThanOrEqualTo")
    O.constraints_thk_A = [O.constraints_thk_A; thk_A];
    O.constraints_thk_b = [O.constraints_thk_b; ...
        constraints.LessThanOrEqualTo];
    return;
end

if isfield(constraints, "GreaterThanOrEqualTo")
    O.constraints_thk_A = [O.constraints_thk_A; -thk_A];
    O.constraints_thk_b = [O.constraints_thk_b; ...
        -constraints.GreaterThanOrEqualTo];
    return;
end

error("One of 'IsFixed', 'LessThanOrEqualTo', or " + ...
    "'GreaterThanOrEqualTo', must be specified.");

end

