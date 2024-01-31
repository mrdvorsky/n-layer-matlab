function [AhHat, AeHat] = computeAhat(O, kRhoP)
%COMPUTEAHAT Computes the matrices AhHat(kRho) and AeHat(kRho).
%   This function computes the matrices AhHat(kRho) and AeHat(kRho) as a
%   function of kRho. The outputs of this function can be used to compute
%   the matrices Ah and Ae, respectively by integrating over kRhoP over the
%   interval [0, 1]. Note that the change of variables from kRhoP to kRho
%   will be applied in this function automatically.
%
% Example Usage:
%   [AhHat, AeHat] = O.computeAhat(kRhoP);
%
% Inputs:
%   kRhoP - A vector of kRhoP coordinates (coordinate transform of kRho).
%       Values should be in the interval [0, 1].
% Outputs:
%   AhHat, AeHat - Matrices computed as a function of kRhoP that can be
%       used to compute the matrix A (i.e., Ah + Ae). The size will be
%       numel(kRhoP) by O.numModes by O.numModes.
%
% Author: Matt Dvorsky

arguments
    O;
    kRhoP(:, 1);
end

%% Integral Change of Variables
% The integral needs to be evaluated from kRho = [0, inf). However, a change
% of variables kRho = L(1 - kRhoP)/kRhoP is used here so that the interpolant
% can be uniform in (0, 1].
L = O.integralScaleFactor;
kRho = L * (1 - kRhoP) ./ kRhoP;

% Weighting function to account for change of variables.
weights_kRho = L ./ (kRhoP.^2);

%% Compute Waveguide Mode Cutoffs
aj(1, 1, :, 1) = O.modesTE(:, 1) * pi ./ O.waveguideA;
bj(1, 1, :, 1) = O.modesTE(:, 2) * pi ./ O.waveguideB;
ai(1, :, 1, 1) = O.modesTE(:, 1) * pi ./ O.waveguideA;
bi(1, :, 1, 1) = O.modesTE(:, 2) * pi ./ O.waveguideB;

%% Helper Functions
f = @(kx, ky) (O.waveguideA.^2 .* O.waveguideB.^2 / pi.^2) ...
    .* sinc((0.5/pi)*O.waveguideA * (kx - aj)) ...
    .* sinc((0.5/pi)*O.waveguideA * (kx - ai)) ...
    .* sinc((0.5/pi)*O.waveguideB * (ky - bj)) ...
    .* sinc((0.5/pi)*O.waveguideB * (ky - bi)) ...
    ./ ((aj + kx) .* (ai + kx) ...
    .* (bj + ky) .* (bi + ky));

%% Compute Integrals Over kPhi at All Values of kRho
% Use 4th dimension for integration over kPhi
[kPhi(1, 1, 1, :), weights_kPhi(1, 1, 1, :)] = ...
    O.fejer2(O.integralPoints_kPhi, 0, 0.5*pi);
kx = kRho .* cos(kPhi);
ky = kRho .* sin(kPhi);

fsc = sum(weights_kPhi .* kRho.^3 .* sin(kPhi).^2 .* cos(kPhi).^2 .* f(kx, ky), 4);
fss = sum(weights_kPhi .* kRho.^3 .* sin(kPhi).^4 .* f(kx, ky), 4);
fcc = sum(weights_kPhi .* kRho.^3 .* cos(kPhi).^4 .* f(kx, ky), 4);

%% Compute AhHat and AeHat Submatrices
AhhhHat = -bj.*bj.*ai.*fsc - aj.*aj.*ai.*fsc;
AhheHat = -aj.*bj.*ai.*fsc + bj.*aj.*ai.*fsc;
AhehHat =  bj.*bj.*bi.*fsc + aj.*aj.*bi.*fsc;
AheeHat =  aj.*bj.*bi.*fsc - bj.*aj.*bi.*fsc;

AehhHat =  bj.*bj.*ai.*fsc - aj.*aj.*ai.*fss;
AeheHat =  aj.*bj.*ai.*fsc + bj.*aj.*ai.*fss;
AeehHat =  bj.*bj.*bi.*fcc - aj.*aj.*bi.*fsc;
AeeeHat =  aj.*bj.*bi.*fcc + bj.*aj.*bi.*fsc;

%% Compute Ah(kRhoP) and Ae(kRhoP)
% Get indices of valid TE and TM modes.
indTE = find(O.modesTE(:, 1) > 0);
indTM = find(O.modesTE(:, 2) > 0);

% Assemble output matrices.
AhHat = weights_kRho .* cat(2, ...
    cat(3, AhhhHat(:, indTE, indTE), AhheHat(:, indTE, indTM)), ...
    cat(3, AhehHat(:, indTM, indTE), AheeHat(:, indTM, indTM)));
AeHat = weights_kRho .* cat(2, ...
    cat(3, AehhHat(:, indTE, indTE), AeheHat(:, indTE, indTM)), ...
    cat(3, AeehHat(:, indTM, indTE), AeeeHat(:, indTM, indTM)));

%% Fix Nans Caused by Singularities At Endpoints
AhHat(kRhoP == 0, :, :) = 0;
AeHat(kRhoP == 0, :, :) = 0;

AhHat(kRhoP == 1, :, :) = 0;
AeHat(kRhoP == 1, :, :) = 0;

end

