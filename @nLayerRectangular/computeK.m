function [kA, kB] = computeK(O, f)
%COMPUTEK Computes kA and kB.
% Computes the matrices kA and kB, which are used to compute the
% unnormalized mode S-parameter matrix.
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

aj(1, :) = O.modesTE(:, 1) * pi ./ O.waveguideA;
bj(1, :) = O.modesTE(:, 2) * pi ./ O.waveguideB;

kjj = conj(sqrt(k0.^2 - aj.^2 - bj.^2));

%% Compute kA and kB Submatrices
% Note that instead of right multiplying by the diagonal matrices KA and
% KB, we can instead multiply each column by the corresponding diagonal
% element for better efficiency. Thus, kA and kB will be constructed as row
% vectors instead of diagonal matrices.

kAhh = k0 + 0*kjj;
kAee = 0*k0 + kjj;

kBhh = 0*k0 + kjj;
kBee = k0 + 0*kjj;

%% Compute kA and kB
% Get indices of valid TE and TM modes.
indTE = find(O.modesTE(:, 1) > 0);
indTM = find(O.modesTE(:, 2) > 0);

% Assemble output vectors. See above note about row vectors vs matrices.
kA = [kAhh(1, indTE, :), kAee(1, indTM, :)];
kB = [kBhh(1, indTE, :), kBee(1, indTM, :)];

end

