function [AhHat, AeHat] = computeAhat(O, kRhoP, modeStruct)
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
    modeStruct;
end

%% Compute Mode Power Levels
[kPhiPower(:, 1), weights_kPhiPower(:, 1)] = ...
    O.fejer2(400, 0, 0.5*pi);
weights_kPhiPower = 4*weights_kPhiPower;

power_TE = ones(numel(O.modeSpectrumTE_Ex), 1);
for ii = 1:numel(power_TE)
    crossProdFun = @(kx, ky, kRho, kPhi) ...
          O.modeSpectrumTE_Ex{ii}(kx, ky, kRho, kPhi).^2 ...
        + O.modeSpectrumTE_Ey{ii}(kx, ky, kRho, kPhi).^2;
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
% power_TM2


end

%% Function to Integrate Over kPhi
function [integralResult] = integralOverKPhi(fun, kRho, kPhi, weights_kPhi)

kx = kRho .* cos(kPhi);
ky = kRho .* sin(kPhi);

% kPhi and weights_kPhi should be vectors along the 2nd dimension.
integralResult = kRho .* sum(fun(kx, ky, kRho, kPhi) .* weights_kPhi, 1);

end

