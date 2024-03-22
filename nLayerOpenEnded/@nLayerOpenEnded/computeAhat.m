function [krcNodes, Ah_weights, Ae_weights] = computeAhat(O)
%COMPUTEAHAT Computes the matrices AhHat(kRho) and AeHat(kRho).
% This function ...
%
% Example Usage:
%   [AhHat, AeHat] = O.computeAhat(kRhoP);
%
% Inputs:
%   kRhoP - A vector of kRhoP coordinates (coordinate transform of kRho).
%       Values should be in the interval [0, 1].
%
% Outputs:
%   AhHat, AeHat - Matrices computed as a function of kRhoP that can be
%       used to compute the matrix A (i.e., Ah + Ae). The size will be
%       numel(kRhoP) by O.numModes by O.numModes.
%
% Author: Matt Dvorsky

arguments
    O;
end

%% Set Scale Factors and Integral Point Counts
% k0MaxAll = [O.modeStructs.MaxOperatingWavenumber];
% k0Max = max(k0MaxAll([O.modeStructs.IsExcitationMode]));
k0Max = max(O.frequencyRange) * 2*pi ./ O.speedOfLight;

L = (4*pi) ./ max([O.modeStructs.ApertureWidth]);
Lc = k0Max;
Lch = Lc;
Lcw = 10*Lc;

Nm = 100;
Nrho = 1*8192;
Nphi = 1*64;

% Nm = O.integral_pointsKrc;
% Nrho = O.integral_pointsKr;
% Nphi = O.integral_points;


% Dimension Assignment:
%   1: moment index
%   2: m
%   3: n
%   4: kr
%   5: kphi

%% Get kRho Weights and Nodes for Integration
% Use 4th dimension for integration over kr
[krcNodes(:, 1), krc, momentH_weights, momentE_weights] = ...
    nLayer.getContourWeights(Nm, Nrho, L, Lc, Lch, Lcw);

%% Compute Weights and Nodes for Integral Over kPhi
% Use 5th dimension for integration over kphi
[kphi(1, 1, 1, 1, :), weights_kphi(1, 1, 1, 1, :)] = ...
    fejer2(Nphi, 0, 0.5*pi);
weights_kphi = 4*weights_kphi;

% if strcmp(modeStruct.ModeSymmetryX, "None")
%     kphi = cat(5, kphi, kphi + 0.5*pi);
%     weights_kphi = 0.5 * cat(5, weights_kphi, weights_kphi);
% end
% 
% if strcmp(modeStruct.ModeSymmetryY, "None")
%     kphi = cat(5, kphi, -flip(kphi));
%     weights_kphi = 0.5 * cat(5, weights_kphi, flip(weights_kphi));
% end

% kphi = cat(5, kphi, kphi + 0.5*pi);
% weights_kphi = 0.5 * cat(5, weights_kphi, weights_kphi);
% 
% kphi = cat(5, kphi, -flip(kphi));
% weights_kphi = 0.5 * cat(5, weights_kphi, flip(weights_kphi));

% if ~any(strcmp([O.modeStructs.SymmetryAxial], "None"))
%     kphi = 0;
%     weights_kphi = 2*pi;
% end

kx = krc .* cos(kphi);
ky = krc .* sin(kphi);

%% Compute Mode Spectrums
modeSpecExm = zeros(1, numel(O.modeStructs), 1, numel(krc), numel(kphi));
modeSpecEym = zeros(1, numel(O.modeStructs), 1, numel(krc), numel(kphi));
for ii = 1:numel(O.modeStructs)
    modeSpecExm(1, ii, 1, :, :) = zeros(size(kx)) ...
        + O.modeStructs(ii).ExSpec(kx, ky, krc, kphi);
    modeSpecEym(1, ii, 1, :, :) = zeros(size(kx)) ...
        + O.modeStructs(ii).EySpec(kx, ky, krc, kphi);
end

offsetX(1, :) = [O.modeStructs.OffsetX];
offsetY(1, :) = [O.modeStructs.OffsetY];
modeSpecExm = modeSpecExm .* exp(-1j .* offsetX .* kx) .* exp(-1j .* offsetY .* ky);
modeSpecEym = modeSpecEym .* exp(-1j .* offsetX .* kx) .* exp(-1j .* offsetY .* ky);

modeSpecExn = circshift(reshape(modeSpecExm, size(modeSpecExm, [1, 3, 2, 4, 5])), 2*Nphi, 5);
modeSpecEyn = circshift(reshape(modeSpecEym, size(modeSpecEym, [1, 3, 2, 4, 5])), 2*Nphi, 5);

%% Compute Cmn
Exy = eye(size(modeSpecExm, 2));

%% Compute Moments
cosPhi = cos(kphi);
sinPhi = sin(kphi);

Ah_moments = innerProduct(momentE_weights, ...
    innerProduct(weights_kphi .* ...
                (sinPhi.*modeSpecExm - cosPhi.*modeSpecEym), ...
                (sinPhi.*modeSpecExn - cosPhi.*modeSpecEyn), ...
                5) .* krc, 4);

Ae_moments = innerProduct(momentH_weights, ...
    innerProduct(weights_kphi .* ...
                (cosPhi.*modeSpecExm + sinPhi.*modeSpecEym), ...
                (cosPhi.*modeSpecExn + sinPhi.*modeSpecEyn), ...
                5) .* krc, 4);

%% Compute Nodes and Weights
Ah_weights = zeros(size(Ah_moments));
Ae_weights = zeros(size(Ae_moments));
for ii = 1:size(Ah_weights(:, :), 2)
    [~, Ah_weights(:, ii)] = fejer2_halfOpen(Nm, Lc, ...
        WeightingMoments=Ah_moments(:, ii));

    [~, Ae_weights(:, ii)] = fejer2_halfOpen(Nm, Lc, ...
        WeightingMoments=Ae_moments(:, ii));
end
Ah_weights = Ah_weights ./ (1 + (krcNodes).^1);
Ae_weights = Ae_weights .* (1 + (krcNodes).^1);

end


