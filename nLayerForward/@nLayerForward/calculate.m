function [gam] = calculate(O, f, er, ur, thk, options)
%CALCULATE Calculate reflection/transmission coefficient(s) for structure.
% Computes the reflection/transmission coefficients when looking into a
% multilayer structure defined by er, ur, and thk at the frequencies
% defined by "f". Some options may depend on the specific subclass of
% nLayerForward (e.g., nLayerRectangular, etc.). Check documentation of the
% specific subclass for more specific information.
%
% Note that the units for all distance and time (or frequency) parameters
% are defined by the speedOfLight parameter. The default units are mm and
% GHz.
%
% Example Usage:
%   NL = nLayerRectangular(maxM, maxN, waveguideBand=wgBand);
%   NL = nLayerCircularTE(numModes, waveguideBand=wgBand);
%   gam = NL.calculate(f, er, ur, thk, ...);
%   gam = NL.calculate(f, {er1, er2}, {}, {thk1, thk2}, ...);
%   gam = NL.calculate(f, {}, ur, thk, ...);
%   gam = NL.calculate(f, er, {}, thk, ..., BackingConductivity=sigma);
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
%   gam - Array of calculated reflection/transmission coefficients for the
%       structure(s). Might be a cell array for each output channel,
%       depending on the specific implementation.
%
% Named Arguments:
%   BackingConductivity (inf) - Conductivity of the backing conductor, in
%       (S/unitDistance). Must have compatible size with all other inputs.
%       A negative conductivity value will be magnetic conductivity.
%
% Author: Matt Dvorsky

arguments
    O;
    f {mustBeNonempty};
    er(:, 1);
    ur(:, 1);
    thk(:, 1) {mustBeNonempty};
    options.BackingConductivity = inf;
end

%% Check Values and Sizes of f, er, ur, and thk
[er, ur, thk] = nLayer.validateStructure(er, ur, thk, ...
    CheckStructureValues=O.checkStructureValues);

%% Change Backing Conductivity
if any(isfinite(options.BackingConductivity(:)))
    [er, ur, thk] = O.changeStructureConductivity(f, er, ur, thk, ...
        options.BackingConductivity);
end

%% Calculate Reflection/Transmission Coefficients
gam = O.calculate_impl(f, er, ur, thk);

end

