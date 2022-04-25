function [er, ur, thk] = changeStructureConductivity(f, er, ur, thk, sigma)
%CALCULATE Modify multilayer structure with new backing conductivity.
% Modifies the structure so that it is conductor-backed with conductivity
% equal to sigma.
%
% Inputs:
%   f - Vector of frequencies (GHz).
%   er - Array of complex relative permittivities for each layer. Every row
%       er(ff, :) should contain the permittivity of each layer at the
%       frequency f(ff). Pass in [] to use default value (1). Optionally,
%       er can be a single row vector.
%   ur - Same as er, except for complex relative permeability.
%   thk - Row vector of thicknesses for each layer. The length of thk
%       should be the same as the number of columns in er and ur. If the
%       last element is inf, this layer will be replaced with an equivalent
%       half-space conductor layer with conductivity sigma. Otherwise, a
%       new layer will be added.
%   sigma - Conductivity value to use, in S/m. Can be a column vector with
%       the same length as f.
% Outputs:
%   er - Array of complex relative permittivities for each layer. Every row
%       er(ff, :) will contain the permittivity of each layer at the
%       frequency f(ff). The last layer value will be set such that it is
%       equivalent to an infinite half-space conductor with conductivity
%       sigma.
%   ur - Same as er, except for complex relative permeability.
%   thk - Row vector of thicknesses for each layer. The length of thk
%       will be the same as the number of columns in er and ur. The last
%       layer will be inf.
%
% Author: Matt Dvorsky

arguments
    f(:, 1);
    er(:, :);
    ur(:, :);
    thk(1, :);
    sigma(:, 1);
end

%% Check Values and Sizes of f, er, ur, and thk
[f, er, ur, thk] = nLayerForward.verifyStructure(f, er, ur, thk, ...
    CheckStructureValues=false);

%% Calculate Equivalent erPrime for Given sigma
erp = sigma ./ ((2*pi .* 8.854e-3) * f);

%% Add or Change Last Layer
if isfinite(thk(end))
    er = [er, 1 - 1j*erp];
    ur = [ur, ones(size(ur, 1), 1)];
    thk = [thk, inf];
else
    er(:, end) = 1 - 1j*erp;
    ur(:, end) = 1;
end

