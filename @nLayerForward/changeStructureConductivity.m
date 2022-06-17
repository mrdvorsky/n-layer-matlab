function [er, ur, thk] = changeStructureConductivity(O, f, er, ur, thk, sigma)
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
%   sigma - Conductivity value to use, in S/unitDistance. Can be a column
%       vector with the same length as f. The distance unit is the same as
%       the units of the speed of light.
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
    O;
    f(:, 1);
    er(:, :);
    ur(:, :);
    thk(1, :);
    sigma(:, 1);
end

%% Check Values and Sizes of f, er, ur, and thk
[er, ur, thk] = nLayerForward.validateStructure(f, er, ur, thk, ...
    CheckStructureValues=false);

%% Calculate Equivalent Loss Factor for Given sigma
mu0_times_c_HenrysPerSecond = 119.9169832 * pi;
k0 = 2*pi .* f ./ O.speedOfLight;

erpp = mu0_times_c_HenrysPerSecond .* sigma ./ k0;

%% Add or Change Last Layer
if isfinite(thk(end))
    er = [er, 1 - 1j*erpp];
    ur = [ur, ones(size(ur, 1), 1)];
    thk = [thk, inf];
else
    er(:, end) = 1 - 1j*erpp;
    ur(:, end) = 1;
end

