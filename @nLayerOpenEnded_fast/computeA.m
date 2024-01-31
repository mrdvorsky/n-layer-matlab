function [A] = computeA(O, f, er, ur, thk)
%COMPUTEA Compute the matrix A for each frequency.
% This function computes the matrix A as a function of each frequency
% specified by f, which is used to compute the unnormalized mode
% S-parameter matrix.
%
% Example Usage (for single frequency, unnormalized S-parameter matrix):
%   [A] = O.computeA(f(ii), er, ur, thk);
%   [B] = O.computeB();
%   [kA, kB] = O.computeK(f(ii));
%   S = (A.*kA + B.*kB) \ (-A.*kA + B.*kB);
%
% Although the example above shows usage with a scalar value for "f", the
% input "f" can be a vector. In this case, the size of the 3rd dimension
% of each output matrix will be equal to numel(f).
%
% Inputs
%   f - Vector of frequencies in GHz.
%   er - Array of complex relative permittivities for each layer (must
%       have compatible dimensions with f).
%   ur - Array of complex relative permeabilities for each layer (must
%       have compatible dimensions with f).
%   thk - Array of thicknesses for each layer (must have compatible 
%       dimensions with er and ur). The last element of thk should have a 
%       value of inf in the infinite halfspace case.
% Outputs:
%   A - Array of O.numModes by O.numModes matrices for each frequency. The
%       size of A will be O.numModes by O.numModes by numel(f).
%
% Author: Matt Dvorsky

arguments
    O;
    f;
    er;
    ur;
    thk;
end

%% Calculate Freespace Wavenumber
k0(1, 1, 1, :) = 2*pi .* f(:) ./ O.speedOfLight;

%% Calculate Gamma0h, Gamma0e
% computeGamma0 expects er and ur to be 1 by numLayers by 1 by length(k0).
[Gamma0h, Gamma0e] = nLayer.computeGamma0(O.fixed_kRho, k0, ...
    permute(er, [3, 2, 4, 1]), permute(ur, [3, 2, 4, 1]), thk);

%% Compute A
% A = innerProduct(Gamma0h, O.fixed_Ah_weights, 1) ...
%     + innerProduct(Gamma0e, O.fixed_Ae_weights, 1);
A = pagemtimes(Gamma0h, "transpose", O.fixed_Ah_weights, "none") ...
    + pagemtimes(Gamma0e, "transpose", O.fixed_Ae_weights, "none");

%% Format Output
A = reshape(A, O.numModes, O.numModes, length(k0));

end

