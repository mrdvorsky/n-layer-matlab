function [A] = computeA(O, f, er, ur, thk)
%COMPUTEA Compute the matrix A for each frequency.
% This function computes the matrix A as a function of each frequency
% specified by f, which is used to compute the unnormalized mode
% S-parameter matrix.
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
% Inputs
%   f - Vector of frequencies in GHz.
%   er - Array of complex relative permittivities for each layer (must
%       have compatible dimensions with f).
%   ur - Array of complex relative permeabilities for each layer (must
%       have compatible dimensions with f).
%   thk - Array of thicknesses for each layer (must have compatible 
%       dimensions with er and ur). The last element of thk should have a 
%       value of inf in the infinite halfspace case.
% Outputs:
%   A - Array of O.numModes by O.numModes matrices for each frequency. The
%       size of A will be O.numModes by O.numModes by numel(f).
%
% Author: Matt Dvorsky

arguments
    O;
    f;
    er;
    ur;
    thk;
end

%% Calculate Freespace Wavenumber
k0(1, 1, 1, :) = 2*pi .* f(:) ./ O.speedOfLight;

%% Initialize A and Gamma0h, Gamma0e
A = zeros(1, O.numModes, O.numModes, length(k0));

% computeGamma0 expects er and ur to be 1 by numLayers by 1 by length(k0).
[Gamma0h, Gamma0e] = O.computeGamma0(O.fixed_kRho, k0, ...
    permute(er, [3, 2, 4, 1]), permute(ur, [3, 2, 4, 1]), thk);

%% Check Error in A When Using Fixed-Point Integration
% Use precomputed weights and nodes for fixed point integration to check
% the accuracy of the dominant mode coefficient. The fixed point weights
% and nodes (i.e., fixed_*) are computed in the "recomputeInterpolants"
% member function.
A_DomMode = sum(Gamma0h .* O.fixed_AhHat(:, 1, 1) ...
    + Gamma0e .* O.fixed_AeHat(:, 1, 1), 1);
errorA_DomMode = sum(Gamma0h .* O.fixed_errorAhHat(:, 1, 1) ...
    + Gamma0e .* O.fixed_errorAeHat(:, 1, 1), 1);

% Determine which frequencies can utilize the fixed point algorithm.
fixedPointFrequency = (abs(errorA_DomMode) ./ abs(A_DomMode)) < O.convergenceAbsTol;

%% Compute A using Fixed Point Integration
% Use the fixed point integration method to compute the integral at all
% frequencies at which this method is sufficiently accurate.
A(:, :, :, fixedPointFrequency) = ...
    sum(Gamma0h(:, 1, 1, fixedPointFrequency) .* O.fixed_AhHat ...
    + Gamma0e(:, 1, 1, fixedPointFrequency) .* O.fixed_AeHat, 1);

% Print message if verbosity flag is 1 or higher.
if O.verbosity > 0
    numLossyFreq = sum(fixedPointFrequency);
    fprintf("Fixed-point integral method used for (%d/%d) frequencies. ", ...
        numLossyFreq, length(k0));
    if numLossyFreq < length(k0)
        fprintf("Using adaptive integral method for the remaining (%d) frequencies.\n", ...
            length(k0) - numLossyFreq);
    else
        fprintf("\n");
    end
end

%% Compute A Using Adaptive Integration
% Use adaptive integration to compute A at the remaining frequencies.
for ff = 1:length(k0)
    if ~fixedPointFrequency(ff)
        % Print message for each frequency if verbosity flag 2 or higher.
        if O.verbosity > 1
            fprintf("f = %.4g GHz - ", f(ff));
        end
        
        A(:, :, :, ff) = O.integralVectorized(...
            @(kRhoP) O.integrandAhat(kRhoP, k0(ff), er(ff, :), ur(ff, :), thk), ...
            0, 1, RelTol=O.convergenceAbsTol, Verbosity=(O.verbosity > 1), ...
            InitialIntervalCount=O.integralInitialSegmentCount);
    end
end

%% Format Output
A = reshape(A, O.numModes, O.numModes, length(k0));

end

