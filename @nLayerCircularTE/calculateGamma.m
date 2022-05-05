function [gam] = calculateGamma(O, f, er, ur, thk)
%CALCULATE Calculate S11 for circular waveguide TE01 mode excitation.
% Computes the reflection coefficient of the circular waveguide TE01
% mode when looking into a multilayer structure defined by er, ur, thk at
% the frequencies defined by f.
%
% Example Usage:
%   NL = nLayerCircularTE(numModes, R=wg_r);
%   gam = NL.calculate(f, er, ur, thk);
%   gam = NL.calculate(f, er, [], thk);
%   gam = NL.calculate(f, [], ur, thk);
%   gam = NL.calculate(f, er, [], thk, BackingConductivity=sigma);
%
% Inputs:
%   f - Column vector of frequencies (GHz).
%   er - Array of complex relative permittivities for each layer. Every row
%       er(ff, :) should contain the permittivity of each layer at the
%       frequency f(ff). Pass in [] to use default value (1).
%   ur - Same as er, except for complex relative permeability.
%   thk - Vector of thicknesses for each layer. The length of thk should be
%       the same as the number of columns in er and ur. Last element can be
%       inf to represent an infinite half-space.
% Outputs:
%   gam - Column vector of reflection coefficients for the TE01 mode. Same
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

%% Get A2
% A2 is precomputed in the "recomputeInterpolants" function.
A2 = O.A2;

%% Assemble Frequency Info (k_A2, k_b2)
[k_A2, k_b2] = O.constructFrequencyMultipliers(f);

%% Calculate Reflection Coefficient at Each Frequency
gam = zeros(length(f), 1);
for ff = 1:length(f)
    x = (A1(:, :, ff) + ur(ff, 1).*A2(:, :).*k_A2(:, :, ff)) ...
        \ (-A1(:, 1, ff) + ur(ff, 1).*A2(:, 1).*k_b2(:, :, ff));
    gam(ff) = x(1);
end

end

