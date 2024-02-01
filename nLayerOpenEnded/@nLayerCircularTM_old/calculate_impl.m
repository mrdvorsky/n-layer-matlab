function [gam] = calculate_impl(O, f, er, ur, thk)
%CALCULATE_IMPL Calculate S11 for circular waveguide TM01 mode excitation.
% Computes the reflection coefficient of the circular waveguide TM01 mode
% when looking into a multilayer structure defined by er, ur, thk at the
% frequencies defined by f.
%
% Inputs:
%   f - Column vector of frequencies.
%   er - Array of complex relative permittivities for each layer. Every row
%       er(ff, :) contains the permittivity of each layer at the frequency
%       f(ff).
%   ur - Same as er, except for complex relative permeability.
%   thk - Vector of thicknesses for each layer. The length of thk is the
%       same as the number of columns in er and ur. Last element can be
%       inf to represent an infinite half-space.
% Outputs:
%   gam - Column vector of reflection coefficients for the TM01 mode. Same
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

%% Convert Cell Array to Vector
er = 0*f + [er{:}];
ur = 0*f + [ur{:}];
thk = [thk{:}];

%% Check for Zero Thickness
if all(thk == 0)
    gam = complex(-ones(size(f)));
    return;
end

%% Compute A, B, kA, and kB
% The call to "O.computeA(...)" is computationally intensive.
[A] = O.computeA(f, er, ur, thk);
[B] = O.computeB();
[kA, kB] = O.computeK(f);

%% Calculate Reflection Coefficient at each Frequency
gam = zeros(length(f), O.numModes, O.numModes);
for ff = 1:length(f)
    Sij = ( A(:, :, ff).*kA(:, :, ff) + B(:, :).*kB(:, :, ff)) ...
        \ (-A(:, :, ff).*kA(:, :, ff) + B(:, :).*kB(:, :, ff));
    gam(ff, :, :) = Sij;
end

gam = gam(:, 1);

end

