function [K] = computeK(O, f)
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

betaCutoff_TE = conj(sqrt(k0.^2 .* O.waveguideEr .* O.waveguideUr - O.cutoffBeta_TE.^2));
betaCutoff_TM = conj(sqrt(k0.^2 .* O.waveguideEr .* O.waveguideUr - O.cutoffBeta_TM.^2));

betaCutoff_TE = complex(real(betaCutoff_TE), -abs(imag(betaCutoff_TE)));
betaCutoff_TM = complex(real(betaCutoff_TM), -abs(imag(betaCutoff_TM)));

%% Compute K
% Assemble output vectors. See above note about row vectors vs matrices.
betaC_TE_over_k0 = betaCutoff_TE ./ (O.waveguideUr .* k0);
betaC_TM_over_k0 = betaCutoff_TM ./ (O.waveguideEr .* k0);

K = sqrt([1./betaC_TE_over_k0; betaC_TM_over_k0]);
K = K .* pagetranspose(K);

end

