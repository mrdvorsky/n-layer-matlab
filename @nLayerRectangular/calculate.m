function [gam] = calculate(O, f, er, ur, thk, options)
%CALCULATE Calculate S11 for TE10 mode excitation.
% Computes the reflection coefficient of the rectangular waveguide TE10
% mode when looking into a multilayer structure defined by er, ur, thk at
% the frequencies defined by f.
%
% Example Usage:
%   NL = nLayerRectangular(maxM, maxN, band=wgBand);
%   gam = NL.calculate(f, er, ur, thk);
%   gam = NL.calculate(f, er, [], thk, AbsTol=1e-4);
%   gam = NL.calculate(f, [], ur, thk);
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
%   gam - Column vector of reflection coefficients for the TE10 mode. Same
%       size as f.
% Named Options:
%   AbsTol - Absolute tolerance used to calculate gam. The default value of
%       this is the member variable convergenceAbsTol.
%
% Author: Matt Dvorsky

arguments
    O;
    f(:, 1);
    er(:, :);
    ur(:, :);
    thk(1, :);
    options.AbsTol(1, 1) = O.convergenceAbsTol;
end

%% Check Values and Sizes of f, er, ur, and thk
[f, er, ur, thk] = O.verifyStructure(f, er, ur, thk, ...
    CheckStructureValues=O.checkStructureValues);

%% Compute A1 and b1
% This is the computationally intensive part of this algorithm
[A1, b1] = O.computeA1b1(f, er, ur, thk, options.AbsTol);

%% Get A2, b2, and etaR1
% A2 and b2 are precomputed in the "recomputeInterpolants" function.
A2 = O.A2;
b2 = O.b2;

etaR1 = sqrt(ur(:, 1) ./ er(:, 1));

%% Assemble Frequency Info (k_A1, k_A2, k_b1, k_b2)
[k_A1, k_A2, k_b1, k_b2] = O.constructFrequencyMultipliers(f);

%% Calculate Reflection Coefficient at each Frequency
gam = zeros(length(f), 1);
for ff = 1:length(f)
    x = (A1(:, :, ff).*k_A1(:, :, ff) + etaR1(ff).*A2.*k_A2(:, :, ff)) ...
        \ (b1(:, :, ff).*k_b1(:, :, ff) + etaR1(ff).*b2.*k_b2(:, :, ff));
    gam(ff) = x(1);
end

end

