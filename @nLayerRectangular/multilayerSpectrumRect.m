function [specE, specH] = multilayerSpectrumRect(tau, k0, er, ur, thk)
%MULTILAYERSPECTRUMRECT Calculate reflection coefficient spectrum.
% This function computes the spectrum for the multilayer structure
% reflection coefficient for a rectangular waveguide. Specifically, it
% computes k_1/(zeta_1*D1^(e)) and zeta_1*D1^(h)/k_1 as a function of tau.
%
% Inputs
%   tau - Spectral wavenumber (can be an array of any size)
%   k0 - Free-space wavenumber (must have compatible dimensions with tau).
%   er - Array of complex relative permittivities for each layer (must have
%       compatible dimensions with k0).
%   ur - Array of complex relative permeabilities for each layer (must have
%       compatible dimensions with k0).
%   thk - Array of thicknesses for each layer (must have compatible 
%       dimensions with er and ur). The last element of thk should have a 
%       value of inf in the infinite halfspace case.
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
k = k0 .* sqrt(er(end) .* ur(end));
zetaPrev = sqrt(k.^2 - tau.^2);
zetaPrev = complex(real(zetaPrev), -abs(imag(zetaPrev)));

if isfinite(thk(end))
    Ce = exp(-2j .* zetaPrev .* thk(end));
else
    Ce = 0;
end
Ch = -Ce;

%% Calculate Structure Reflection Coefficient
% Loop over layers starting with second to last
for ii = (length(thk) - 1):-1:1
    k = k0 .* sqrt(er(ii) .* ur(ii));
    zeta = sqrt(k.^2 - tau.^2);
    zeta = complex(real(zeta), -abs(imag(zeta)));
    
    Be = er(ii + 1) ./ er(ii) .* zeta ./ zetaPrev;
    Bh = ur(ii + 1) ./ ur(ii) .* zeta ./ zetaPrev;
    
    Ce = exp(-2j .* zeta .* thk(ii)) .* (Ce .* (1 + Be) - (1 - Be)) ...
        ./ (-Ce .* (1 - Be) + (1 + Be));
    Ch = exp(-2j .* zeta .* thk(ii)) .* (Ch .* (1 + Bh) - (1 - Bh)) ...
        ./ (-Ch .* (1 - Bh) + (1 + Bh));
    
    zetaPrev = zeta;
end

%% Calculate Output
specE = k .* (1 + Ce) ./ (1 - Ce) ./ zetaPrev;
specH = zetaPrev .* (1 - Ch) ./ (1 + Ch) ./ k;

end

