function [Smn] = calculate(O, f, er, ur, thk)
%CALCULATE Calculate Smn for a multimoded waveguide looking into structure.
% Computes the S-parameter matrix (Smn) of an open-ended multimoded
% waveguide looking into a multilayer structure defined by "er", "ur",
% "thk", and at the frequencies defined by "f".
%
% Example Usage:
%   NL = nLayerRectangular(maxM, maxN, waveguideBand="Ka");
%   NL = nLayerCircular(0, numModes, waveguideBand="Ka_TE01", ...
%           modeSymmetryAxial="TE");
%
%   gam = NL.calculate(f, er, ur, thk);
%   gam = NL.calculate(f, {er1, er2}, {}, {thk1, thk2});
%   gam = NL.calculate(f, {}, ur, thk);
%
% Inputs:
%   f - Array of frequencies. Must have compatible size with each layer of
%       "er" and "ur", but this is not checked.
%   er - Cell array of complex relative permittivities for each layer.
%       Every element of the cell array corresponds to one layer of the
%       structure, and each must be a compatible size to "f". For example,
%       "er{n}" corresponds to the nth layer. Pass in {} to use the default
%       value of 1.
%   ur - Same as "er", except for complex relative permeability.
%   thk - Same as "er" and "ur" but for the thicknesses of each layer.
%       Obviously, the value of "thk" should not change with frequency, but
%       this is not checked.
%
% Outputs:
%   Smn - The computed mode S-parameter matrix.
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

