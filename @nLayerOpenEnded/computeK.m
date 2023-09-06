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

betaC_TE = conj(sqrt(k0.^2 .* O.waveguideEr .* O.waveguideUr - O.modeBetaCutoffTE.^2));
betaC_TM = conj(sqrt(k0.^2 .* O.waveguideEr .* O.waveguideUr - O.modeBetaCutoffTM.^2));

betaC_TE = complex(real(betaC_TE), -abs(imag(betaC_TE)));
betaC_TM = complex(real(betaC_TM), -abs(imag(betaC_TM)));

%% Compute kA and kB Submatrices
% Note that instead of right multiplying by the diagonal matrices KA and
% KB, we can instead multiply each column by the corresponding diagonal
% element for better efficiency. Thus, kA and kB will be constructed as row
% vectors instead of diagonal matrices.



%% Compute kA and kB
% Assemble output vectors. See above note about row vectors vs matrices.
betaC_TE_over_k0 = betaC_TE ./ (O.waveguideUr .* k0);
betaC_TM_over_k0 = betaC_TM ./ k0;


kA = sqrt([1./betaC_TE_over_k0, betaC_TM_over_k0]);
kA = kA .* pagetranspose(kA);
kB = kA*0 + eye(O.numModes);

end

