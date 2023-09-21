function [Smn] = calculate_impl(O, f, er, ur, thk)
%CALCULATE_IMPL Calculate S11 for rectangular waveguide TE10 mode excitation.
% Computes the reflection coefficient of the rectangular waveguide TE10
% mode when looking into a multilayer structure defined by er, ur, thk at
% the frequencies defined by f.
%
% Inputs:
%   f - Column vector of frequencies.
%   er - Array of complex relative permittivities for each layer. Every row
%       er(ff, :) contains the permittivity of each layer at the frequency
%       f(ff).
%   ur - Same as er, except for complex relative permeability.
%   thk - Vector of thicknesses for each layer. The length of thk is the
%       same as the number of columns in er and ur. Last element can be
%       inf to represent an infinite half-space.
% Outputs:
%   gam - Column vector of reflection coefficients for the TE10 mode. Same
%       size as f.
%
% Author: Matt Dvorsky

arguments
    O;
    f(:, 1);
    er(:, :);
    ur(:, :);
    thk(1, :);
end

%% Check for Zero Thickness
if all(thk == 0)
    Smn = complex(-ones(size(f)));
    return;
end

%% Compute A, B, kA, and kB
% The call to "O.computeA(...)" is computationally intensive.
[A] = O.computeA(f, er, ur, thk);
[K] = O.computeK(f);

%% Calculate S-parameter Matrix at each Frequency
A_times_K = A.*K;
idMat = eye(size(A, 1));

Smn = pagemldivide(idMat + A_times_K, idMat - A_times_K);

%% Get S-parameter Submatrix
outputModes_All = [O.outputModes_TE; O.outputModes_TM; O.outputModes_Hybrid];
Smn = permute(Smn(outputModes_All, outputModes_All, :), [3, 1, 2]);

end

