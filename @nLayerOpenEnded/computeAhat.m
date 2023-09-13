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

%% Integral Change of Variables
% The integral needs to be evaluated from kRho = [0, inf). However, a change
% of variables kRho = L(1 - kRhoP)/kRhoP is used here so that the interpolant
% can be uniform in (0, 1].
L = modeStruct.IntegralScaleFactor;
kRho = L * (1 - kRhoP) ./ kRhoP;

% Fix kRho to never be infinite or zero. The integrands at these endpoints
% will be set to zero later, so the specific values don't matter.
kRho(kRhoP == 0) = nan;
kRho(kRhoP == 1) = nan;
kRho(kRhoP == 0) = max(kRho, [], "omitnan") * 1.01;
kRho(kRhoP == 1) = min(kRho, [], "omitnan") * 0.99;

% Weighting function to account for change of variables.
weights_kRho = L ./ (kRhoP.^2);

%% Compute Weights and Nodes for Integral Over kPhi
% Use 4th dimension for integration over kPhi
[kPhi(1, 1, 1, :), weights_kPhi(1, 1, 1, :)] = ...
    O.fejer2(O.integralPoints_kPhi, 0, 0.5*pi);
weights_kPhi = 4*weights_kPhi;

if strcmp(modeStruct.ModeSymmetryX, "None")
    kPhi = cat(4, kPhi, kPhi + 0.5*pi);
    weights_kPhi = 0.5 * cat(4, weights_kPhi, weights_kPhi);
end

if strcmp(modeStruct.ModeSymmetryY, "None")
    kPhi = cat(4, kPhi, -flip(kPhi));
    weights_kPhi = 0.5 * cat(4, weights_kPhi, flip(weights_kPhi));
end

if ~strcmp(modeStruct.ModeSymmetryAxial, "None")
    kPhi = 0;
    weights_kPhi = 2*pi;
end

kx = kRho .* cos(kPhi);
ky = kRho .* sin(kPhi);

%% Compute Mode Spectrums
% Calculate TE mode spectrums.
modeSpecEx_TE = zeros(numel(kRho), numel(modeStruct.SpecEx_TE), 1, numel(kPhi));
modeSpecEy_TE = zeros(numel(kRho), numel(modeStruct.SpecEy_TE), 1, numel(kPhi));
for ii = 1:numel(modeStruct.SpecEx_TE)
    modeSpecEx_TE(:, ii, 1, :) = zeros(size(kx)) ...
        + modeStruct.SpecEx_TE{ii}(kx, ky, kRho, kPhi);
    modeSpecEy_TE(:, ii, 1, :) = zeros(size(kx)) ...
        + modeStruct.SpecEy_TE{ii}(kx, ky, kRho, kPhi);
end

% Calculate TM mode spectrums.
modeSpecEx_TM = zeros(numel(kRho), numel(modeStruct.SpecEx_TM), 1, numel(kPhi));
modeSpecEy_TM = zeros(numel(kRho), numel(modeStruct.SpecEy_TM), 1, numel(kPhi));
for ii = 1:numel(modeStruct.SpecEx_TM)
    modeSpecEx_TM(:, ii, 1, :) = zeros(size(kx)) ...
        + modeStruct.SpecEx_TM{ii}(kx, ky, kRho, kPhi);
    modeSpecEy_TM(:, ii, 1, :) = zeros(size(kx)) ...
        + modeStruct.SpecEy_TM{ii}(kx, ky, kRho, kPhi);
end

% Calculate Hybrid mode spectrums.
modeSpecEx_Hybrid = zeros(numel(kRho), numel(modeStruct.SpecEx_Hybrid), 1, numel(kPhi));
modeSpecEy_Hybrid = zeros(numel(kRho), numel(modeStruct.SpecEy_Hybrid), 1, numel(kPhi));
for ii = 1:numel(modeStruct.SpecEx_Hybrid)
    modeSpecEx_Hybrid(:, ii, 1, :) = zeros(size(kx)) ...
        + modeStruct.SpecEx_Hybrid{ii}(kx, ky, kRho, kPhi);
    modeSpecEy_Hybrid(:, ii, 1, :) = zeros(size(kx)) ...
        + modeStruct.SpecEy_Hybrid{ii}(kx, ky, kRho, kPhi);
end

%% Combine Modes
modeSpecExm = cat(2, modeSpecEx_TE, modeSpecEx_TM, modeSpecEx_Hybrid);
modeSpecEym = cat(2, modeSpecEy_TE, modeSpecEy_TM, modeSpecEy_Hybrid);

modeSpecExn = reshape(modeSpecExm, size(modeSpecExm, [1, 3, 2, 4]));
modeSpecEyn = reshape(modeSpecEym, size(modeSpecEym, [1, 3, 2, 4]));

%% Compute Ahat
AhHat = weights_kRho .* innerProduct(weights_kPhi .* ...
    (ky.*modeSpecExm - kx.*modeSpecEym), ...
    (ky.*modeSpecExn - kx.*modeSpecEyn), ...
    4) ./ kRho;

AeHat = weights_kRho .* innerProduct(weights_kPhi .* ...
    (kx.*modeSpecExm + ky.*modeSpecEym), ...
    (kx.*modeSpecExn + ky.*modeSpecEyn), ...
    4) ./ kRho;

%% Fix Nans Caused by Singularities At Endpoints
AhHat(kRhoP == 0, :, :) = 0;
AeHat(kRhoP == 0, :, :) = 0;

AhHat(kRhoP == 1, :, :) = 0;
AeHat(kRhoP == 1, :, :) = 0;

end


