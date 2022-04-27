function [specH] = multilayerSpectrumH(tau, k0, er, ur, thk)
%MULTILAYERSPECTRUMH Calculate reflection coefficient spectrum.
% This function computes the spectrum for the multilayer structure
% reflection coefficient for a rectangular waveguide. Specifically, it
% computes k_1/(zeta_1*D1^(e)) and zeta_1*D1^(h)/k_1 as a function of tau.
%
% Inputs
%   tau - Column vector of spectral wavenumbers.
%   k0 - Array of free-space wavenumbers (the size must be 1-by-1-by-...,
%       where the 3rd and 4th dimensions may be any size).
%   er - Array of complex relative permittivities for each layer (the
%       size must be 1-by-numLayers-by-..., where the 3rd and 4th
%       dimensions must be compatible with k0.
%   ur - Same as er, except for complex relative permeabilities.
%   thk - Vector of thicknesses for each layer (must have same length as
%       size(er, 2) and size(ur, 2). The last element of thk should have a
%       value of inf for the infinite halfspace case.
% Outputs
%   specH - Calculated spectrum, same size as (tau .* k0)
%
% The output of this function can be used along with the output of
% "computeIntegrandH" to compute the integral I(m, n). See documentation
% for "computeIntegrandH" for more details.
%
% The output specH is equal to zeta_1*D1^(h) as a function of tau.
%
% Author: Matt Dvorsky

%% Calculate Last Layer
zetaPrev = sqrt(k0.^2 .* er(1, end, :, :) .* ur(1, end, :, :) - tau.^2);
zetaPrev = complex(real(zetaPrev), -abs(imag(zetaPrev)));

if isfinite(thk(end))
    Ch = -exp(-2j .* zetaPrev .* thk(end));
else
    Ch = 0;
end

%% Calculate Structure Reflection Coefficient
% Loop over layers starting with second to last
for ii = (length(thk) - 1):-1:1
    zeta = sqrt(k0.^2 .* er(1, ii, :, :) .* ur(1, ii, :, :) - tau.^2);
    zeta = complex(real(zeta), -abs(imag(zeta)));
    
    Bh = ur(1, ii + 1, :, :) ./ ur(1, ii, :, :) .* zeta ./ zetaPrev;
    
    Ch = exp(-2j .* zeta .* thk(ii)) .* (Ch .* (1 + Bh) - (1 - Bh)) ...
        ./ (-Ch .* (1 - Bh) + (1 + Bh));
    
    zetaPrev = zeta;
end

%% Calculate Output
specH = zetaPrev .* (1 - Ch) ./ (1 + Ch);

end

