function [er, ur, thk] = validateStructure(er, ur, thk, options)
%Check the multilayer structure and correct any issues.
% This function can be used to verify proper dimensions and values of the
% inputs er, ur, thk.
%
% Inputs:
%   er - Cell array of complex relative permittivities for each layer.
%       Every element of the cell array corresponds to one layer of the
%       structure, and each must be compatible sizes. For example, "er{n}"
%       corresponds to the nth layer. Pass in {} to use the default value
%       of 1. If it is a vector, will be converted to cell array.
%   ur - Same as "er", except for complex relative permeability.
%   thk - Same as "er" and "ur".
%
% Outputs:
%   er - Same as input "er", but will be a cell array of arrays.
%   ur - Same as input "ur", but will be a cell array of arrays.
%   thk - Same as input "thk", but will be a cell array of arrays.
%
% Named Arguments:
%   CheckStructureValues (true) - If true, check that "er", "ur", and
%       "thk" contain physical values.
%
% Author: Matt Dvorsky

arguments
    er(:, 1);
    ur(:, 1);
    thk(:, 1) {mustBeNonempty};

    options.CheckStructureValues(1, 1) logical = true;
    options.RequireConstantValuesPerLayer(1, 1) logical = false;
end

%% Check for Vector Inputs
if ~iscell(er)
    er = mat2cell(er, ones(numel(er), 1));
end

if ~iscell(ur)
    ur = mat2cell(ur, ones(numel(ur), 1));
end

if ~iscell(thk)
    thk = mat2cell(thk, ones(numel(thk), 1));
end

%% Check for Empty or Singleton er and ur
if isempty(er)
    er = repmat({1}, numel(thk), 1);
end

if isempty(ur)
    ur = repmat({1}, numel(thk), 1);
end

%% Check Lengths
if numel(er) ~= numel(ur) || numel(er) ~= numel(thk)
    error("Inputs 'er', 'ur', and 'thk' must all be cell " + ...
        "arrays with the same length (or empty).");
end

%% Check that Each Layer is Scalar
if options.RequireConstantValuesPerLayer ...
        && (any(cellfun(@numel, er) ~= 1) ...
        ||  any(cellfun(@numel, ur) ~= 1) ...
        ||  any(cellfun(@numel, thk) ~= 1))
    error("The inputs 'er', 'ur', and 'thk' must be " + ...
        "constant for each layer.");
end

%% Check Structure Values ("er", "ur", "thk")
if ~options.CheckStructureValues
    return;
end

% Check value finiteness.
for n = 1:numel(thk) - 1
    if ~all(isfinite(er{n}(:))) || ~all(isfinite(ur{n}(:))) ...
            || ~all(isfinite(thk{n}(:)))
        error("All elements of 'er', 'ur', and 'thk' except the " + ...
            "last layer must be finite.");
    end
end

% Check passivity.
for n = 1:numel(thk)
    if ~all(real(er{n}(:)) >= 1, "all") || ~all(real(ur{n}(:)) >= 1, "all")
        error("The real parts of 'er' and 'ur' must be greater than 1. " + ...
            "To disable this check, set 'checkStructureValues' to false.");
    end
    if ~all(imag(er{n}(:)) <= 0, "all") || ~all(imag(ur{n}(:)) <= 0, "all")
        error("The imaginary parts of 'er' and 'ur' must be nonpositive. " + ...
            "To disable this check, set 'checkStructureValues' to false.");
    end
    if ~all(thk{n}(:) >= 0)
        error("All elements of 'thk' must be nonnegative. " + ...
            "To disable this check, set 'checkStructureValues' to false.");
    end
end

end

