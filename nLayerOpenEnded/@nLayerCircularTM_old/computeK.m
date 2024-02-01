function [kA, kB] = computeK(O, f)
%COMPUTEK Computes kA and kB.
% Computes the matrices kA and kB, which are used to compute the
% normalized mode S-parameter matrix.
%
% Inputs:
%   f - vector of frequencies to consider.
% Outputs:
%   kA, kB - Matrices used to compute S-parameters.
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
% Author: Matt Dvorsky

arguments
    O;
    f(1, 1, :) double;
end

%% Mode Coefficients
k0 = 2*pi .* f ./ O.speedOfLight;

kc0j(1, :) = O.modeCutoffs;

k0j = conj(sqrt(k0.^2 - kc0j.^2));

%% Compute kA and kB
% Note that instead of right multiplying by the diagonal matrices KA and
% KB, we can instead multiply each column by the corresponding diagonal
% element for better efficiency. Thus, kA and kB will be constructed as row
% vectors instead of diagonal matrices.

kA = 0*k0 + k0j;
kB = k0 + 0*k0j;

end

