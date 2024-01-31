function [gam] = calculate(O, f, er, ur, thk, varargin, options)
%CALCULATE Calculate reflection/trasnmission coefficient(s) for structure.
% Computes the reflection/transmission coefficients when looking into a
% multilayer structure defined by er, ur, and thk at the frequencies
% defined by f. Some options may depend on the specific subclass of
% nLayerForward (e.g., nLayerRectangular, etc.). Check documentation of the
% specific subclass for more specific information.
%
% Note that the units for all distance and time (or frequency) parameters
% are defined by the speedOfLight parameter. The default units are mm and
% GHz.
%
% Example Usage:
%   NL = nLayerRectangular(maxM, maxN, waveguideBand=wgBand);
%   NL = nLayerCircularTE(numModes, waveguideR=wgR);
%   gam = NL.calculate(f, er, ur, thk, ...);
%   gam = NL.calculate(f, er, [], thk, ...);
%   gam = NL.calculate(f, [], ur, thk, ...);
%   gam = NL.calculate(f, [], ur, thk, ..., BackingConductivity=sigma);
%
% Inputs:
%   f - Column vector of frequencies.
%   er - Array of complex relative permittivities for each layer. Every row
%       er(ff, :) should contain the permittivity of each layer at the
%       frequency f(ff). Pass in [] to use default value (1).
%   ur - Same as er, except for complex relative permeability.
%   thk - Vector of thicknesses for each layer. The length of thk should be
%       the same as the number of columns in er and ur. Last element can be
%       inf to represent an infinite half-space.
% Outputs:
%   gam - Column vector(s) of reflection coefficients for the TE10 mode.
%       First dimension will have the same size as f.
% Named Arguments:
%   BackingConductivity (inf) - Conductivity of the backing conductor, in
%       (S/unitDistance). Must be scalar or have the same length as f.
%
% Author: Matt Dvorsky

arguments
    O;
    f(:, 1);
    er(:, :);
    ur(:, :);
    thk(1, :);
end
arguments (Repeating)
    varargin;
end
arguments
    options.BackingConductivity(:, 1) = inf;
end

%% Check Values and Sizes of f, er, ur, and thk
[er, ur, thk] = O.validateStructure(f, er, ur, thk, ...
    CheckStructureValues=O.checkStructureValues);

%% Change Backing Conductivity
if any(isfinite(options.BackingConductivity))
    [er, ur, thk] = O.changeStructureConductivity(f, er, ur, thk, ...
        options.BackingConductivity);
end

%% Calculate Reflection/Transmission Coefficients
gam = O.calculate_impl(f, er, ur, thk, varargin{:});

end

