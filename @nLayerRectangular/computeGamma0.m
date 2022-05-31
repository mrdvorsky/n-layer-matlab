function [GammaH, GammaE] = computeGamma0(kRho, k0, er, ur, thk)
%COMPUTEGAMMA0 Calculate reflection coefficient spectrum.
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
%   specE - Calculated spectrum, same size as (tau .* k0)
%   specH - Calculated spectrum, same size as (tau .* k0)
%
% The output of this function can be used along with the output of
% "computeIntegrandEH" to compute the integrals I^(e)_ii(m, n, p, q) and
% I^(h)_ii(m, n, p, q). See documentation for "computeIntegrandEH" for
% more details.
%
% The outputs specE and specH are equal to k_1/(zeta_1*D1^(e)) and
% zeta_1*D1^(h)/k_1, respectively, as a function of tau.
%
% Author: Matt Dvorsky

%% Calculate Last Layer
kzPlus1 = sqrt(k0.^2 .* er(1, end, :, :) .* ur(1, end, :, :) - kRho.^2);
kzPlus1 = complex(real(kzPlus1), -abs(imag(kzPlus1)));

if isfinite(thk(end))
    GammaE = exp(-2j .* kzPlus1 .* thk(end));
else
    GammaE = 0;
end
GammaH = -GammaE;

%% Calculate Gamma_i for Each Layer
% Loop over layers starting with second to last
for ii = (length(thk) - 1):-1:1
    kz = sqrt(k0.^2 .* er(1, ii, :, :) .* ur(1, ii, :, :) - kRho.^2);
    kz = complex(real(kz), -abs(imag(kz)));
    
    gammaH = (ur(1, ii + 1, :, :).*kz - ur(1, ii, :, :).*kzPlus1) ...
        ./ (ur(1, ii + 1, :, :).*kz + ur(1, ii, :, :).*kzPlus1);
    gammaE = (er(1, ii + 1, :, :).*kz - er(1, ii, :, :).*kzPlus1) ...
        ./ (er(1, ii + 1, :, :).*kz + er(1, ii, :, :).*kzPlus1);
    
    GammaH = exp(-2j .* kz .* thk(ii)) .* (gammaH + GammaH) ...
        ./ (1 + gammaH.*GammaH);
    GammaE = exp(-2j .* kz .* thk(ii)) .* (gammaE + GammaE) ...
        ./ (1 + gammaE.*GammaE);
    
    kzPlus1 = kz;
end

%% Calculate Gamma_0
GammaE = k0 .* er(1, 1, :, :) .* (1 + GammaE) ./ ((1 - GammaE) .* kzPlus1);
GammaH = kzPlus1 .* (1 - GammaH) ./ ((1 + GammaH) .* k0 .* ur(1, 1, :, :));

end

