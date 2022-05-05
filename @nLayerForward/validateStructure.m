function [er, ur, thk] = validateStructure(f, er, ur, thk, options)
%VALIDATESTRUCTURE Check the multilayer structure and correct any issues.
% This function can be used to verify proper dimensions and values of the
% inputs f, er, ur, thk.
%
% Inputs:
%   f - Vector of frequencies (GHz).
%   er - Array of complex relative permittivities for each layer. Every row
%       er(ff, :) should contain the permittivity of each layer at the
%       frequency f(ff). Pass in [] to use default value (1). Optionally,
%       er can be a single row vector or a scalar value.
%   ur - Same as er, except for complex relative permeability.
%   thk - Row vector of thicknesses for each layer. The length of thk
%       should be the same as the number of columns in er and ur. Last
%       element can be inf to represent an infinite half-space.
% Outputs:
%   er - Array of complex relative permittivities for each layer. Every row
%       er(ff, :) will contain the permittivity of each layer at the
%       frequency f(ff).
%   ur - Same as er, except for complex relative permeability.
%   thk - Row vector of thicknesses for each layer. The length of thk
%       will be the same as the number of columns in er and ur.
%
% Author: Matt Dvorsky

arguments
    f(:, 1) {mustBeNonempty};
    er(:, :);
    ur(:, :);
    thk(1, :) {mustBeNonempty};
    options.CheckStructureValues {mustBeNumericOrLogical} = true;
end

%% Check for Empty or Singleton er and ur
if isempty(er)
    er = ones(length(f), length(thk));
elseif numel(er) == 1
    er = er + zeros(size(thk));
end

if isempty(ur)
    ur = ones(length(f), length(thk));
elseif numel(ur) == 1
    ur = ur + zeros(size(thk));
end

%% Check Values of er, ur, and thk
if options.CheckStructureValues
    if ~all(real(er) >= 1, "all") || ~all(real(ur) >= 1, "all")
        error("The real parts of er and ur must be greater than 1.");
    end
    if ~all(imag(er) <= 0, "all") || ~all(imag(ur) <= 0, "all")
        error("The imaginary parts of er and ur must be nonpositive.");
    end
    if ~all(thk >= 0)
        error("All elements of thk must be nonnegative.");
    end
    if ~all(isfinite(thk(1:end - 1)))
        error("All elements of thk except the last must be finite.");
    end
end

%% Check Dimensions of f, er, ur, and thk
if ~all(size(er, [1, 2]) == [length(f), length(thk)])
    if isvector(er) && (size(er, 2) == length(thk))
        er = repmat(er, length(f), 1);
    else
        error(strcat("The size of er should be either 1-by-numLayers ", ...
            "or numFreqs-by-numLayers."));
    end
end

if ~all(size(ur, [1, 2]) == [length(f), length(thk)])
    if isvector(ur) && (size(ur, 2) == length(thk))
        ur = repmat(ur, length(f), 1);
    else
        error(strcat("The size of ur should be either 1-by-numLayers ", ...
            "or numFreqs-by-numLayers."));
    end
end

