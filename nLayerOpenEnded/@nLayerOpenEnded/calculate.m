function [Smn] = calculate(O, f, er, ur, thk)
%CALCULATE Calculate Smn for a multimoded waveguide looking into structure.
% Computes the S-parameter matrix (Smn) of an open-ended multimoded
% waveguide looking into a multilayer structure defined by "er", "ur",
% "thk", and at the frequencies defined by "f".
%
% Inputs:
%
% Outputs:
%
% Author: Matt Dvorsky

arguments
    O;
    f(:, 1);
    er(1, :);
    ur(1, :);
    thk(1, :) {mustBeNonempty};
end

%% Check if Integral Weights Should be Regenerated
if O.shouldRecomputeWeights
    O.computeIntegralWeights();
end

%% Validate Structure
[er, ur, thk] = nLayer.validateStructure(er, ur, thk, ...
    CheckStructureValues=O.checkStructureValues);

%% Compute A and K
[A] = O.computeA(f, er, ur, thk);
[K] = O.computeK(f);

%% Calculate S-parameter Matrix at each Frequency
A_times_K = A.*K;
idMat = eye(size(A, 1));

Smn = pagemldivide(idMat + A_times_K, idMat - A_times_K);

%% Get S-parameter Submatrix
Smn = permute(Smn(O.receiveModeIndices, O.excitationModeIndices, :), [3, 1, 2]);

end

