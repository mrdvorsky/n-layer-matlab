function [ specE, specH ] = multilayerSpectrumRect( tau, k0, er, ur, thk )
%MULTILAYERSPECTRUMTE Calculate multilayer reflection coefficient spectrum
%   for circular waveguide TM01 nLayer
% Inputs
%   tau - spectral wavenumber (can be an array of any size)
%   k0 - free-space wavenumber (must have compatible dimensions with tau)
%   er - vector of complex relative permittivities for each layer
%   ur - vector of complex relative permeabilities for each layer
%   thk - vector of thicknesses for each layer (same length as er and ur)
% Outputs
%   specE - Calculated spectrum, same size as (tau .* k0)
%   specH - Calculated spectrum, same size as (tau .* k0)

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

