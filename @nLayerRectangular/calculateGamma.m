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
A1 = O.computeA1(f, er, ur, thk);

%% Get A2 and etaR1
% A2 is precomputed in the "recomputeInterpolants" function.
A2 = O.A2;

etaR1 = sqrt(ur(:, 1) ./ er(:, 1));

%% Assemble Frequency Info (k_A1, k_A2, k_b1, k_b2)
[k_A1, k_A2, k_b1, k_b2] = O.constructFrequencyMultipliers(f);

%% Calculate Reflection Coefficient at each Frequency
gam = zeros(length(f), 1);
for ff = 1:length(f)
    x = (A1(:, :, ff).*k_A1(:, :, ff) + etaR1(ff).*A2(:, :).*k_A2(:, :, ff)) ...
        \ (-A1(:, 1, ff).*k_b1(:, :, ff) + etaR1(ff).*A2(:, 1).*k_b2(:, :, ff));
    gam(ff) = x(1);
end

end

