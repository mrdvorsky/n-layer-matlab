function [Gamma0h, Gamma0e, Gamma0hp, Gamma0ep] = computeGammaA(kr, k0, er, ur, thk)
%Computes Gamma0h and Gamma0e for the multilayer structure.
% This function computes the spectrum for the multilayer structure
% reflection coefficient for TE and TM modes as a function of "kr".
% Specifically, it computes "Gamma0h(kr)" and "Gamma0e(kr)".
%
% Inputs
%   kr - Array of target spectral wavenumbers.
%   k0 - Array of free-space wavenumbers. Must have compatible size
%       with "kr" (and all other inputs).
%   er - Cell array of complex relative permittivities for each layer.
%       Every element of the cell array corresponds to one layer of the
%       structure, and each must be a compatible size to "k0". For example,
%       "er{n}" corresponds to the nth layer. Pass in {} to use the default
%       value of 1.
%   ur - Same as "er", except for complex relative permeability.
%   thk - Same as "er" and "ur" but for the thicknesses of each layer.
%       Obviously, the value of "thk" should not change with frequency, but
%       this is not checked. The last layer can be infinite.
%
% Outputs
%   Gamma0h - Calculated spectrum, the size will be the compatible size of
%       all inputs.
%   Gamma0e - Calculated spectrum, same size as "Gamma0h".
%
% Author: Matt Dvorsky

%% Calculate Last Layer
kzPlus1 = sqrt(k0.^2 .* er{end} .* ur{end} - kr.^2);
kzPlus1 = complex(real(kzPlus1), -abs(imag(kzPlus1)));

if all(~isfinite(thk{end}(:)))
    GammaE = 0;
else
    GammaE = exp(-2j .* kzPlus1 .* thk{end});
    GammaE(isnan(GammaE)) = 0;
end
GammaH = -GammaE;

%% Calculate Gamma_i for Each Layer
% Loop over layers starting with second to last
for ii = (numel(thk) - 1):-1:1
    kz = sqrt(k0.^2 .* er{ii} .* ur{ii} - kr.^2);
    kz = complex(real(kz), -abs(imag(kz)));
    
    gammaH = (ur{ii + 1} ./ ur{ii}) .* (kz ./ kzPlus1);
    gammaE = (er{ii + 1} ./ er{ii}) .* (kz ./ kzPlus1);

    gammaH(isnan(gammaH)) = 0;
    gammaE(isnan(gammaE)) = 0;
    
    GammaH = exp(-2j .* kz .* thk{ii}) ...
        .* (GammaH .* (gammaH + 1) + (gammaH - 1)) ...
        ./ (GammaH .* (gammaH - 1) + (gammaH + 1));
    GammaE = exp(-2j .* kz .* thk{ii}) ...
        .* (GammaE .* (gammaE + 1) + (gammaE - 1)) ...
        ./ (GammaE .* (gammaE - 1) + (gammaE + 1));
    
    kzPlus1 = kz;
end

%% Calculate Gamma0h and Gamma0e
% Gamma0h = kzPlus1 .* (1 - GammaH) ./ ((1 + GammaH) .* k0 .* ur{1});
% Gamma0e = k0 .* er{1} .* (1 + GammaE) ./ ((1 - GammaE) .* kzPlus1);

Gamma0h = GammaH + zeros(size(kzPlus1));
Gamma0e = GammaE + zeros(size(kzPlus1));

Gamma0hp = kzPlus1 .* (1 - GammaH) ./ (k0 .* ur{1});
Gamma0ep = k0 .* er{1} .* (1 + GammaE) ./ kzPlus1;

end

