function [krcNodes, Ah_weights, Ae_weights] = computeAhat(O, modeStruct)
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
    modeStruct;
end

%% Set Scale Factors and Integral Point Counts
L = 5;
Lch = 4;
Lcw = 40;

Nm = 300;
Nrho = 8192;
Nphi = 64;

% Dimension Assignment:
%   1: moment index
%   2: m
%   3: n
%   4: kr
%   5: kphi

%% Get kRho Weights and Nodes for Integration
% Use 4th dimension for integration over kr
[krcNodes(:, 1), krc, momentH_weights, momentE_weights] = ...
    nLayer.getContourWeights(Nm, Nrho, L, Lch, Lcw);

%% Compute Weights and Nodes for Integral Over kPhi
% Use 5th dimension for integration over kphi
[kphi(1, 1, 1, 1, :), weights_kphi(1, 1, 1, 1, :)] = ...
    fejer2(Nphi, 0, 0.5*pi);
weights_kphi = 4*weights_kphi;

if strcmp(modeStruct.ModeSymmetryX, "None")
    kphi = cat(4, kphi, kphi + 0.5*pi);
    weights_kphi = 0.5 * cat(4, weights_kphi, weights_kphi);
end

if strcmp(modeStruct.ModeSymmetryY, "None")
    kphi = cat(4, kphi, -flip(kphi));
    weights_kphi = 0.5 * cat(4, weights_kphi, flip(weights_kphi));
end

if ~strcmp(modeStruct.ModeSymmetryAxial, "None")
    kphi = 0;
    weights_kphi = 2*pi;
end

kx = krc .* cos(kphi);
ky = krc .* sin(kphi);

%% Compute Mode Spectrums
% Calculate TE mode spectrums.
modeSpecEx_TE = zeros(1, numel(modeStruct.SpecEx_TE), 1, numel(krc), numel(kphi));
modeSpecEy_TE = zeros(1, numel(modeStruct.SpecEx_TE), 1, numel(krc), numel(kphi));
for ii = 1:numel(modeStruct.SpecEx_TE)
    modeSpecEx_TE(1, ii, 1, :, :) = zeros(size(kx)) ...
        + modeStruct.SpecEx_TE{ii}(kx, ky, krc, kphi);
    modeSpecEy_TE(1, ii, 1, :, :) = zeros(size(kx)) ...
        + modeStruct.SpecEy_TE{ii}(kx, ky, krc, kphi);
end

% Calculate TM mode spectrums.
modeSpecEx_TM = zeros(1, numel(modeStruct.SpecEx_TM), 1, numel(krc), numel(kphi));
modeSpecEy_TM = zeros(1, numel(modeStruct.SpecEx_TM), 1, numel(krc), numel(kphi));
for ii = 1:numel(modeStruct.SpecEx_TM)
    modeSpecEx_TM(1, ii, 1, :, :) = zeros(size(kx)) ...
        + modeStruct.SpecEx_TM{ii}(kx, ky, krc, kphi);
    modeSpecEy_TM(1, ii, 1, :, :) = zeros(size(kx)) ...
        + modeStruct.SpecEy_TM{ii}(kx, ky, krc, kphi);
end

% Calculate Hybrid mode spectrums.
modeSpecEx_Hybrid = zeros(1, numel(modeStruct.SpecEx_Hybrid), 1, numel(krc), numel(kphi));
modeSpecEy_Hybrid = zeros(1, numel(modeStruct.SpecEx_Hybrid), 1, numel(krc), numel(kphi));
for ii = 1:numel(modeStruct.SpecEx_Hybrid)
    modeSpecEx_Hybrid(1, ii, 1, :, :) = zeros(size(kx)) ...
        + modeStruct.SpecEx_Hybrid{ii}(kx, ky, krc, kphi);
    modeSpecEy_Hybrid(1, ii, 1, :, :) = zeros(size(kx)) ...
        + modeStruct.SpecEy_Hybrid{ii}(kx, ky, krc, kphi);
end

%% Combine Modes
modeSpecExm = cat(2, modeSpecEx_TE, modeSpecEx_TM, modeSpecEx_Hybrid);
modeSpecEym = cat(2, modeSpecEy_TE, modeSpecEy_TM, modeSpecEy_Hybrid);

modeSpecExn = reshape(modeSpecExm, size(modeSpecExm, [1, 3, 2, 4, 5]));
modeSpecEyn = reshape(modeSpecEym, size(modeSpecEym, [1, 3, 2, 4, 5]));

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
    [~, Ah_weights(:, ii)] = fejer2_halfOpen(Nm, L, ...
        WeightingMoments=Ah_moments(:, ii));

    [~, Ae_weights(:, ii)] = fejer2_halfOpen(Nm, L, ...
        WeightingMoments=Ae_moments(:, ii));
end
Ah_weights = Ah_weights ./ sqrt(1 + (krcNodes).^2);
Ae_weights = Ae_weights .* sqrt(1 + (krcNodes).^2);

%% Check Mode Scaling and Orthogonality
% if modeStruct.CheckModeScalingAndOrthogonality
%     modeCrossMatrix = squeeze(mean(weights_kRho .* sum(weights_kphi ...
%         .* (modeSpecExm.*conj(modeSpecExn) + modeSpecEym.*conj(modeSpecEyn)) ...
%         .* kRho, 4), 1));
% 
%     modeCrossMatrixError = squeeze(mean(weights_kRho .* sum(errorWeights_kPhi ...
%         .* (modeSpecExm.*conj(modeSpecExn) + modeSpecEym.*conj(modeSpecEyn)) ...
%         .* kRho, 4), 1))
% 
%     if any(abs(modeCrossMatrix - eye(size(modeCrossMatrix))) > 0.001, "all")
%         warning("One or more modes may be scaled incorrectly, " + ...
%             "or one or more pairs of modes may not be orthogonal. " + ...
%             "The mode cross inner product matrix is shown below, " + ...
%             "which should be approximately equal to the identity " + ...
%             "matrix. \n\nModeCrossInnerProductMatrix =\n\n" + ...
%             "%s\n Note that the S-parameter matrix will be " + ...
%             "incorrect if this issue is not fixed. To disable " + ...
%             "this warning, set the 'CheckModeScalingAndOrthogonality' " + ...
%             "flag to false in the 'defineWaveguideModes' function.\n", ...
%             formattedDisplayText(modeCrossMatrix));
%     end
% end

end


