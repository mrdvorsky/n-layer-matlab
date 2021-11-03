function [A1, b1] = computeA1b1(O, f, er, ur, thk, AbsTol)
%COMPUTEA1B1 Compute the integral of the matrix A1 and b1 over frequency.
% This function computes the matrices A1 and b1 at each frequency specified
% by f. The output of this function can along with the outputs of the 
% constructFrequencyMultipliers(...) function and the 
% constructMatrixEquation(...) function to calculate S11 for a rectangular
% waveguide using the example code.
%
% Example Usage:
%   [k_A1, k_A2, k_b1, k_b2] = O.constructFrequencyMultipliers(f);
%   [~, A2, ~, b2] = O.constructMatrixEquation();
%   [A1, b1] = computeA1b1(f, er, ur, thk, AbsTol);
%   sqrt(ur(1) ./ er(1));
%   x = (A1.*k_A1 + etaR1.*A2.*k_A2) ...
%        \ (b1.*k_b1 + etaR1.*b2.*k_b2);
%   S11 = x(1);
%
% Inputs
%   f - Vector of frequencies in GHz.
%   er - Array of complex relative permittivities for each layer (must have
%       compatible dimensions with f).
%   ur - Array of complex relative permeabilities for each layer (must have
%       compatible dimensions with k0).
%   thk - Array of thicknesses for each layer (must have compatible 
%       dimensions with er and ur). The last element of thk should have a 
%       value of inf in the infinite halfspace case.
%   AbsTol - Tolerance with which to compute the matrix A1.
% Outputs:
%   A1 - Array of O.numModes by O.numModes matrices for each frequency. The
%       size of A1 will be O.numModes by O.numModes by numel(f).
%   b1 - Array of O.numModes by 1 matrices for each frequency. The size of 
%       b1 will be O.numModes by 1 by numel(f).
%
% Note that the negative of b1 is always equal to the first column of A1.
% This fact will be utilized here.
%
% Author: Matt Dvorsky

arguments
    O;
    f;
    er;
    ur;
    thk;
    AbsTol;
end

%% Calculate Freespace Wavenumber
k0(1, 1, 1, :) = 2*pi .* f(:) ./ O.c;

%% Initialize A1 and specE, specH
A1 = zeros(1, O.numModes, O.numModes, length(k0));

% multilayerSpectrumRect expects er and ur to be
% 1-by-numLayers-by-1-by-length(k0).
specH = O.multilayerSpectrumH(O.fixed_tau, k0, ...
    permute(er, [3, 2, 4, 1]), permute(ur, [3, 2, 4, 1]), thk);

%% Check Error in A1 When Using Fixed-Point Integration
% Use precomputed weights and nodes for fixed point integration to check
% the accuracy of the dominant mode coefficient. The fixed point weights
% and nodes (i.e., fixed_*) are computed in the "recomputeInterpolants"
% member function.
A1DomMode = sum(specH .* O.fixed_A1_H(:, 1, 1), 1);
errorA1DomMode = sum(specH .* O.fixed_errA1_H(:, 1, 1), 1);

% All frequencies that have a low error bound are "lossy".
isLossyFrequency = (abs(errorA1DomMode) ./ abs(A1DomMode)) < AbsTol;

%% Compute A1 at Lossy Frequencies Using Fixed Point Integration
% Use the fixed point integration method to compute the integral at all
% frequencies at which this method is sufficiently accurate (i.e, at
% "lossy" frequencies).
A1(:, :, :, isLossyFrequency) = sum(...
    specH(:, 1, 1, isLossyFrequency) .* O.fixed_A1_H, 1);

% Print message if verbosity flag is 1 or higher
if O.verbosity > 0
    numLossyFreq = sum(isLossyFrequency);
    fprintf("Fixed-point integral method used for (%d/%d) frequencies. ", ...
        numLossyFreq, length(k0));
    if numLossyFreq < length(k0)
        fprintf("Using adaptive integral method for the remaining (%d) frequencies.\n", ...
            length(k0) - numLossyFreq);
    else
        fprintf("\n");
    end
end

%% Compute A1 at Low Loss Frequencies
% Use adaptive integration to compute A1 at the remaining frequencies.
for ff = 1:length(k0)
    if ~isLossyFrequency(ff)
        % Print message for each frequency if verbosity flag 2 or higher
        if O.verbosity > 1
            fprintf("f = %.4g GHz - ", f(ff));
        end
        
        A1(:, :, :, ff) = O.integralVectorized(...
            @(tauP) O.integrandA1(tauP, k0(ff), er(ff, :), ur(ff, :), thk), ...
            0, 1, RelTol=AbsTol, Verbosity=(O.verbosity > 1), ...
            InitialIntervalCount=O.integralInitialSegmentCount);
    end
end

%% Format Outputs
% Note that (-b1) is always the first column of A1.
A1 = reshape(A1, O.numModes, O.numModes, length(k0));
b1 = -A1(:, 1, :);

end

