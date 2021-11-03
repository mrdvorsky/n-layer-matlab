function gam = calculateWithConductivity(O, f, er, ur, thk, sigma, varargin)
%CALCULATEWITHCONDUCTIVITY Calculate gam using a given backing conductivity.
% Computes the reflection and/or transmission coefficient(s) of the when
% looking into a multilayer structure defined by er, ur, thk at the
% frequencies defined by f, after modifying the structure to be backed with
% the specified conductivity values. Uses the calculate function of the
% object used to call this function after modifying the structure.
%
% Example Usage:
%   NL = *nLayerForwardSubclass*(...);
%   gam = NL.calculate(f, er, ur, thk, sigma);
%   gam = NL.calculate(f, er, ur, thk, sigma, Prop1=val1, ...);
%
% Inputs:
%   f - Column vector of frequencies (GHz).
%   er - Array of complex relative permittivities for each layer. Every row
%       er(ff, :) should contain the permittivity of each layer at the
%       frequency f(ff). Pass in [] to use default value (1).
%   ur - Same as er, except for complex relative permeability.
%   thk - Vector of thicknesses for each layer. The length of thk should be
%       the same as the number of columns in er and ur. If the last element
%       is inf, this layer will be replaced with an equivalent half-space
%       conductor layer with conductivity sigma. Otherwise, a new layer
%       will be added.
% Outputs:
%   gam - Column vector(s) of reflection and/or transmission coefficients
%       for the TE10 mode. Each column is the same size as f.
% Named Options: Any named options specified are passed through to the
%   calculate() function of the calling object.
%
% Author: Matt Dvorsky

arguments
    O;
    f(:, 1);
    er(:, :);
    ur(:, :);
    thk(1, :);
    sigma(:, 1)
end

arguments (Repeating)
    varargin;
end

%% Modify Structure with New Conductivity Value
[er, ur, thk] = O.changeStructureConductivity(f, er, ur, thk, sigma);

%% Calculate gam
gam = O.calculate(f, er, ur, thk, varargin{:});

end

