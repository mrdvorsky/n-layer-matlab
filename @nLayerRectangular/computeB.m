function [B] = computeB(O)
%COMPUTEB Computes the matrix B.
%   This function computes the matrix B, which is used to compute
%   the unnormalized mode S-parameter matrix.
%
% Example Usage (for single frequency, unnormalized S-parameter matrix):
%   [A] = O.computeA(f(ii), er, ur, thk);
%   [B] = O.computeB();
%   [kA, kB] = O.computeK(f(ii));
%   S = (A.*kA + B.*kB) \ (-A.*kA + B.*kB);
%
% Outputs:
%   B - The matrix B, which will have size O.numModes by O.numModes.
%
% Author: Matt Dvorsky

arguments
    O;
end

%% Compute Mode Coefficients
ai(:, 1) = O.modesTE(:, 1) * pi ./ O.waveguideA;
bi(:, 1) = O.modesTE(:, 2) * pi ./ O.waveguideB;

%% Compute B Submatrices
Bhh = -diag(ai .* (1 + (bi == 0)));
Bhe =  diag(bi);
Beh =  diag(bi .* (1 + (bi == 0)));
Bee =  diag(ai);

%% Compute B
% Get indices of valid TE and TM modes.
indTE = find(O.modesTE(:, 1) > 0);
indTM = find(O.modesTE(:, 2) > 0);

% Assemble output matrix.
B = (0.25 * O.waveguideA * O.waveguideB) .* [...
    Bhh(indTE, indTE), Bhe(indTE, indTM); ...
    Beh(indTM, indTE), Bee(indTM, indTM)];

end

