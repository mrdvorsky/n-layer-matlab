function [Gamma0h] = computeGamma0(kRho, k0, er, ur, thk)
%COMPUTEGAMMA0 Computes Gamma0h for the multilayer structure.
% This function computes the spectrum for the multilayer structure
% reflection coefficient for TE modes as a function of kRho.
% Specifically, it computes Gamma0h(kRho).
% 
% The function inputs requirements are specifically formatted for nLayer
% waveguide computations.
%
% Inputs
%   kRho - Column vector of spectral wavenumbers.
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
%   Gamma0h - Calculated spectrum, same size as (kRho .* k0)
%
% Author: Matt Dvorsky

%% Calculate Last Layer
kzPlus1 = sqrt(k0.^2 .* er(1, end, :, :) .* ur(1, end, :, :) - kRho.^2);
kzPlus1 = complex(real(kzPlus1), -abs(imag(kzPlus1)));

if isfinite(thk(end))
    GammaH = -exp(-2j .* kzPlus1 .* thk(end));
else
    GammaH = 0;
end

%% Calculate Gamma_i for Each Layer
% Loop over layers starting with second to last
for ii = (length(thk) - 1):-1:1
    kz = sqrt(k0.^2 .* er(1, ii, :, :) .* ur(1, ii, :, :) - kRho.^2);
    kz = complex(real(kz), -abs(imag(kz)));
    
    gammaH = (ur(1, ii + 1, :, :).*kz - ur(1, ii, :, :).*kzPlus1) ...
        ./ (ur(1, ii + 1, :, :).*kz + ur(1, ii, :, :).*kzPlus1);
    
    GammaH = exp(-2j .* kz .* thk(ii)) .* (gammaH + GammaH) ...
        ./ (1 + gammaH.*GammaH);
    
    kzPlus1 = kz;
end

%% Calculate Gamma_0
Gamma0h = kzPlus1 .* (1 - GammaH) ./ ((1 + GammaH) .* ur(1, 1, :, :));

end

