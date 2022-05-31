function [gam] = calculateGamma(O, f, er, ur, thk)
%CALCULATE Calculate S11 for rectangular waveguide TE10 mode excitation.
% Computes the reflection coefficient of the rectangular waveguide TE10
% mode when looking into a multilayer structure defined by er, ur, thk at
% the frequencies defined by f.
%
% Inputs:
%   f - Column vector of frequencies (GHz).
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
    gam = complex(-ones(size(f)));
    return;
end

%% Compute A1
% This is the computationally intensive part of this algorithm
A = O.computeA(f, er, ur, thk);

%% Get A2 and etaR1
% A2 is precomputed in the "recomputeInterpolants" function.
B = O.A2;

etaR1 = sqrt(ur(:, 1) ./ er(:, 1));

%% Assemble Frequency Info (kA, kB)
[kA, kB] = O.computeK(f);

%% Calculate Reflection Coefficient at each Frequency
gam = zeros(length(f), 1);
for ff = 1:length(f)
    S = (A(:, :, ff).*kA(:, :, ff) + etaR1(ff).*B(:, :).*kB(:, :, ff)) ...
        \ (-A(:, 1, ff).*kA(:, :, ff) + etaR1(ff).*B(:, 1).*kB(:, :, ff));
    gam(ff) = S(1);
end

end

