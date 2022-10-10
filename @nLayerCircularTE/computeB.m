function [B] = computeB(O)
%COMPUTEB Computes the matrix B.
%   This function computes the matrix B, which is used to compute
%   the normalized mode S-parameter matrix.
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

%% Compute Mode Cutoffs
kc0i(:, 1) = O.modeCutoffs;

%% Compute B
B = 0.5 * kc0i .* diag(besselj(0, O.waveguideR * kc0i));

end

