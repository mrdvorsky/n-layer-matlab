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

%% Compute Integrals Over kPhi at All Values of kRho
% Use 4th dimension for integration over kPhi
[kPhi(1, 1, 1, :), weights_kPhi(1, 1, 1, :)] = ...
    O.fejer2(O.integralPoints_kPhi, 0, 0.5*pi);
weights_kPhi = 4*weights_kPhi;
kx = kRho .* cos(kPhi);
ky = kRho .* sin(kPhi);

%% Compute Mode Power Levels
[kPhiPower(:, 1), weights_kPhiPower(:, 1)] = ...
    O.fejer2(400, 0, 0.5*pi);
weights_kPhiPower = 4*weights_kPhiPower;

power_TE = ones(numel(O.modeSpectrumTE_Hx), 1);
for ii = 1:numel(power_TE)
    crossProdFun = @(kx, ky, kRho, kPhi) ...
          O.modeSpectrumTE_Hx{ii}(kx, ky, kRho, kPhi).^2 ...
        + O.modeSpectrumTE_Hy{ii}(kx, ky, kRho, kPhi).^2;
    Bfun = @(kRho) integralOverKPhi(crossProdFun, kRho, kPhiPower, weights_kPhiPower);
%     power_TE(ii) = integral(Bfun, 0, inf, RelTol=1e-6, AbsTol=0);
    
%     [nodes(1, :), weights(1, :)] = clenshawCurtisHalfOpen(1024, 1);
    power_TE2(ii) = sum(Bfun(kRho.') .* weights_kRho.', "omitnan") ./ 4096;
end
power_TE2

power_TM = ones(numel(O.modeSpectrumTM_Ex), 1);
for ii = 1:numel(power_TM)
    crossProdFun = @(kx, ky, kRho, kPhi) ...
          O.modeSpectrumTM_Ex{ii}(kx, ky, kRho, kPhi).^2 ...
          + O.modeSpectrumTM_Ey{ii}(kx, ky, kRho, kPhi).^2;
    Bfun = @(kRho) integralOverKPhi(crossProdFun, kRho, kPhiPower, weights_kPhiPower);
%     power_TM(ii) = integral(Bfun, 0, inf, RelTol=1e-6, AbsTol=0);

    power_TM2(ii) = sum(Bfun(kRho.') .* weights_kRho.', "omitnan") ./ 4096;
end
power_TM2

%% Compute Mode Spectrums
% Calculate TE mode Hertzian potentials.
specHxm = zeros(numel(kRho), numel(O.modeSpectrumTE_Hx), 1, numel(kPhi));
specHym = zeros(numel(kRho), numel(O.modeSpectrumTE_Hy), 1, numel(kPhi));
for ii = 1:numel(O.modeSpectrumTE_Hx)
    specHxm(:, ii, 1, :) = O.modeSpectrumTE_Hx{ii}(kx, ky, kRho, kPhi) ...
        + zeros(size(kx));
    specHym(:, ii, 1, :) = O.modeSpectrumTE_Hy{ii}(kx, ky, kRho, kPhi) ...
        + zeros(size(kx));
end
specHxn = reshape(specHxm, size(specHxm, [1, 3, 2, 4]));
specHyn = reshape(specHym, size(specHym, [1, 3, 2, 4]));

% Calculate TM mode Hertzian potentials.
specExm = zeros(numel(kRho), numel(O.modeSpectrumTM_Ex), 1, numel(kPhi));
specEym = zeros(numel(kRho), numel(O.modeSpectrumTM_Ey), 1, numel(kPhi));
for ii = 1:numel(O.modeSpectrumTM_Ex)
    specExm(:, ii, 1, :) = O.modeSpectrumTM_Ex{ii}(kx, ky, kRho, kPhi) ...
        ./ sqrt(abs(power_TM(ii))) + zeros(size(kx));
    specEym(:, ii, 1, :) = O.modeSpectrumTM_Ey{ii}(kx, ky, kRho, kPhi) ...
        ./ sqrt(abs(power_TM(ii))) + zeros(size(kx));
end
specExn = reshape(specExm, size(specExm, [1, 3, 2, 4]));
specEyn = reshape(specEym, size(specEym, [1, 3, 2, 4]));

%% Compute Ahat
AhhhHat = innerProduct(weights_kPhi .* (ky.*specHym + kx.*specHxm), (ky.*specHyn + kx.*specHxn), 4) ./ kRho;
AhheHat = innerProduct(weights_kPhi .* (ky.*specHym + kx.*specHxm), (ky.*specExn - kx.*specEyn), 4) ./ kRho;
AhehHat = innerProduct(weights_kPhi .* (ky.*specExm - kx.*specEym), (ky.*specHyn + kx.*specHxn), 4) ./ kRho;
AheeHat = innerProduct(weights_kPhi .* (ky.*specExm - kx.*specEym), (ky.*specExn - kx.*specEyn), 4) ./ kRho;

AehhHat = innerProduct(weights_kPhi .* (kx.*specHym - ky.*specHxm), (kx.*specHyn - ky.*specHxn), 4) ./ kRho;
AeheHat = innerProduct(weights_kPhi .* (kx.*specHym - ky.*specHxm), (kx.*specExn + ky.*specEyn), 4) ./ kRho;
AeehHat = innerProduct(weights_kPhi .* (kx.*specExm + ky.*specEym), (kx.*specHyn - ky.*specHxn), 4) ./ kRho;
AeeeHat = innerProduct(weights_kPhi .* (kx.*specExm + ky.*specEym), (kx.*specExn + ky.*specEyn), 4) ./ kRho;

%% Compute B
Bhh = mean(sum(weights_kPhi .* specHxm .* specHxn, 4) .* kRho, 1, "omitnan") ...
    + mean(innerProduct(weights_kPhi .* specHym, specHyn, 4) .* kRho, 1, "omitnan");
Bhe = mean(innerProduct(weights_kPhi .* specHym, specExn, 4) .* kRho, 1, "omitnan") ...
    - mean(innerProduct(weights_kPhi .* specHxm, specEyn, 4) .* kRho, 1, "omitnan");
Beh = mean(innerProduct(weights_kPhi .* specHym, specExn, 4) .* kRho, 1, "omitnan") ...
    - mean(innerProduct(weights_kPhi .* specHxm, specEyn, 4) .* kRho, 1, "omitnan");
Bee = mean(innerProduct(weights_kPhi .* specExm, specExn, 4) .* kRho, 1, "omitnan") ...
    + mean(innerProduct(weights_kPhi .* specEym, specEyn, 4) .* kRho, 1, "omitnan");

% O.matrix_B = [squeeze(Bhh), squeeze(Bhe); ...
%         squeeze(Beh), squeeze(Bee)];
% matB = O.matrix_B

%% Assemble Output Matrices
AhHat = weights_kRho .* cat(2, ...
    cat(3, AhhhHat, AhheHat), ...
    cat(3, AhehHat, AheeHat));
AeHat = weights_kRho .* cat(2, ...
    cat(3, AehhHat, AeheHat), ...
    cat(3, AeehHat, AeeeHat));

% AhHat = weights_kRho .* cat(2, ...
%     cat(3, O.waveguideUr.*AhhhHat, AhheHat), ...
%     cat(3, O.waveguideUr.*AhehHat, AheeHat));
% AeHat = weights_kRho .* cat(2, ...
%     cat(3, O.waveguideUr.*AehhHat, AeheHat), ...
%     cat(3, O.waveguideUr.*AeehHat, AeeeHat));

%% Fix Nans Caused by Singularities At Endpoints
AhHat(kRhoP == 0, :, :) = 0;
AeHat(kRhoP == 0, :, :) = 0;

AhHat(kRhoP == 1, :, :) = 0;
AeHat(kRhoP == 1, :, :) = 0;

end

%% Function to Integrate Over kPhi
function [integralResult] = integralOverKPhi(fun, kRho, kPhi, weights_kPhi)

kx = kRho .* cos(kPhi);
ky = kRho .* sin(kPhi);

% kPhi and weights_kPhi should be vectors along the 2nd dimension.
integralResult = kRho .* sum(fun(kx, ky, kRho, kPhi) .* weights_kPhi, 1);

end

