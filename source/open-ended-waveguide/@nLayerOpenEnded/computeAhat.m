function [krcNodes, Ah_weights, Ae_weights] = computeAhat(O)
%Computes the matrices AhHat(kRho) and AeHat(kRho).
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

% Dimension Assignment:
%   1: moment index
%   2: m
%   3: n
%   4: kr
%   5: kphi

%% Set Scale Factors and Integral Point Counts
k0Max = max(O.frequencyRange) * 2*pi ./ O.speedOfLight;

% L = (4*pi) ./ max([O.modeStructs.ApertureWidth]);
Lc = k0Max;
Lch = 0.5*Lc;
Lcw = 10*Lc;
LcFinal = 200*Lc;

krcWaypoints = {0, 0.05*Lch*(1 + 1j), Lch*(1 + 1j), (Lcw + Lch*1j), Lcw + Lch, LcFinal};
a = krcWaypoints(1:end-1);
b = krcWaypoints(2:end);

Nm = O.integral_pointsKrc;
Nrho = O.integral_pointsKr;

Nphi = 4*64;

%% Get kRho Weights and Nodes for Integration
% Use 4th dimension for integration over kr
[krc, moment_weights] = ...
    getContourWeights(Nm, Nrho, a, b);

%% Compute Weights and Nodes for Integral over kphi
% Use 5th dimension for integration over kphi
[kphi(1, 1, 1, 1, :), weights_kphi(1, 1, 1, 1, :)] = ...
    trap(Nphi, 0, 2*pi);
% weights_kphi = 4*weights_kphi;

%% Compute Moment Integrals 1
% for ii = 1:numel(krc)
%     kx = krc{ii} .* cos(kphi);
%     ky = krc{ii} .* sin(kphi);
% 
%     % Compute Mode Spectrums
%     modeSpecExm = zeros(1, numel(O.modeStructs), 1, numel(krc{ii}), numel(kphi));
%     modeSpecEym = zeros(1, numel(O.modeStructs), 1, numel(krc{ii}), numel(kphi));
%     for mm = 1:numel(O.modeStructs)
%         modeSpecExm(1, mm, 1, :, :) = zeros(size(kx)) ...
%             + O.modeStructs(mm).ExSpec(kx, ky, krc{ii}, kphi);
%         modeSpecEym(1, mm, 1, :, :) = zeros(size(kx)) ...
%             + O.modeStructs(mm).EySpec(kx, ky, krc{ii}, kphi);
%     end
% 
%     modeSpecExn = reshape(modeSpecExm, size(modeSpecExm, [1, 3, 2, 4, 5]));
%     modeSpecEyn = reshape(modeSpecEym, size(modeSpecEym, [1, 3, 2, 4, 5]));
% 
%     % Compute Moments
%     cosPhi = cos(kphi);
%     sinPhi = sin(kphi);
% 
%     Ah_moments{ii} = innerProduct(moment_weights{ii}, ...
%         innerProduct(weights_kphi .* ...
%                     (sinPhi.*modeSpecExm - cosPhi.*modeSpecEym), ...
%                     circshift(sinPhi.*modeSpecExn - cosPhi.*modeSpecEyn, 32, 5), ...
%                     5) .* krc{ii}.^(ii>0), 4);
% 
%     Ae_moments{ii} = innerProduct(moment_weights{ii}, ...
%         innerProduct(weights_kphi .* ...
%                     (cosPhi.*modeSpecExm + sinPhi.*modeSpecEym), ...
%                     circshift(cosPhi.*modeSpecExn + sinPhi.*modeSpecEyn, 32, 5), ...
%                     5) .* krc{ii}.^(ii>0), 4);
% end

%% Compute Moment Integrals
% for ii = 1:numel(krc)
%     krc{ii} = single(krc{ii});
% end
% kphi = single(kphi);

for ii = 1:numel(krc)
    kx = krc{ii} .* cos(kphi);
    ky = krc{ii} .* sin(kphi);
    
    % Compute Mode Spectrums
    modeSpecWhm = zeros(1, O.numModes, 1, numel(krc{ii}), numel(kphi));
    modeSpecWem = zeros(1, O.numModes, 1, numel(krc{ii}), numel(kphi));
    for mm = 1:O.numModes
        modeSpecWhm(1, mm, 1, :, :) = zeros(size(kx)) ...
            + O.waveguideModes(mm).WhSpec(kx, ky, krc{ii}, kphi);
        modeSpecWem(1, mm, 1, :, :) = zeros(size(kx)) ...
            + O.waveguideModes(mm).WeSpec(kx, ky, krc{ii}, kphi);
    end
    
    modeSpecWhn = reshape(modeSpecWhm, size(modeSpecWhm, [1, 3, 2, 4, 5]));
    modeSpecWen = reshape(modeSpecWem, size(modeSpecWem, [1, 3, 2, 4, 5]));
    
    Ah_moments{ii} = -innerProduct(moment_weights{ii}, ...
        innerProduct(weights_kphi .* modeSpecWhm, circshift(modeSpecWhn, 128, 5), 5) ...
        .* krc{ii}, 4);
    
    Ae_moments{ii} = -innerProduct(moment_weights{ii}, ...
        innerProduct(weights_kphi .* modeSpecWem, circshift(modeSpecWen, 128, 5), 5) ...
        .* krc{ii}, 4);
end

%% Compute Nodes and Weights
for ii = 1:numel(krc)
    Ah_weights{ii} = zeros(size(Ah_moments{ii}));
    Ae_weights{ii} = zeros(size(Ae_moments{ii}));
    for mm = 1:size(Ah_weights{ii}(:, :), 2)
        [krcNodes{ii}, Ah_weights{ii}(:, mm)] = fejer2(Nm{ii}, a{ii}, b{ii}, ...
            WeightingMoments=Ah_moments{ii}(:, mm));

        [krcNodes{ii}, Ae_weights{ii}(:, mm)] = fejer2(Nm{ii}, a{ii}, b{ii}, ...
            WeightingMoments=Ae_moments{ii}(:, mm));
    end
end

%% Combine Weights and Nodes
krcNodes = cell2mat(krcNodes.');
Ah_weights = cell2mat(Ah_weights.');
Ae_weights = cell2mat(Ae_weights.');

end


