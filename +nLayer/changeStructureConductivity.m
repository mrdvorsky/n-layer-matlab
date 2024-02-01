function [er, ur, thk] = changeStructureConductivity(f, er, ur, thk, sigma, options)
%CHANGESTRUCTURECONDUCTIVITY Modify multilayer structure with new backing conductivity.
% Modifies the structure so that it is conductor-backed with conductivity
% equal to "sigma", in (S/unitDistance). If "sigma" is negative, magnetic
% conductivity will be used instead of electric.
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
%   sigma - Conductivity of the backing conductor, in (S/unitDistance).
%       Must be scalar or have compatible size with "f". A negative
%       conductivity value will be magnetic conductivity. Must have
%       compatible size with all other inputs.
%
% Outputs:
%   er - Same as input "er", except the last layer value will be set such
%       that it is equivalent to an infinite half-space conductor with
%       electric conductivity sigma (if sigma is positive).
%   ur - Same as input "ur", except the last layer value will be set such
%       that it is equivalent to an infinite half-space conductor with
%       magnetic conductivity sigma (if sigma is negative).
%   thk - Same as input "thk", except if the last layer is not "inf", it
%       an additional layer with "inf" will be added.
%
% Named Arguments:
%   SpeedOfLight (299.792458) - Speed of light to used for calculation. All
%       units must be consistent with this value.
%
% Author: Matt Dvorsky

arguments
    f {mustBeNonempty};
    er(:, 1);
    ur(:, 1);
    thk(:, 1) {mustBeNonempty};
    sigma;
    options.SpeedOfLight(1, 1) {mustBePositive} = 299.792458;
end

%% Check Inputs
nLayer.validateStructure(f, er, ur, thk);

%% Calculate Equivalent Loss Factor for Given "sigma"
mu0_times_c_HenrysPerSecond = 119.9169832 * pi;
k0 = 2*pi .* f ./ options.SpeedOfLight;

erpp = mu0_times_c_HenrysPerSecond .* abs(sigma) ./ k0;
urpp = (sign(sigma) < 0) * erpp;
erpp = (sign(sigma) > 0) * erpp;

%% Add or Change Last Layer
if all(isfinite(thk{end}(:)))
    % Add conductor layer.
    thk = [thk; inf];
    er = [er, 1 - 1j*erpp];
    ur = [ur, 1 - 1j*urpp];
elseif all(~isfinite(thk{end}(:)))
    % Conductor layer exists. Modify it.
    er{end} = 1 - 1j*erpp;
    ur{end} = 1 - 1j*urpp;
else
    error("Last layer thickness must be either all finite or all infinite.");
end

