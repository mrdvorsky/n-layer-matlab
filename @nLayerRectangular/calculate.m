function [gam] = calculate(O, f, er, ur, thk, AbsTol)
%CALCULATE Calculate S11 for TE10 mode
% Inputs
%   f - vector of frequencies (GHz)
%   er - vector of complex relative permittivities for each layer
%   ur - vector of complex relative permeabilities for each layer
%   thk - vector of thicknesses for each layer (same length as er and ur)
%   numPoints - number of points to use in the integral quadrature
% Outputs
%   gam - vector of TE10 reflection coefficients over frequency

%% Check Inputs
arguments
    O
    f(:, 1)
    er(1, :)
    ur(1, :)
    thk(1, :)
    AbsTol(1, 1) = O.convergenceAbsTol;
end

% Make sure er, ur, and thk are in the correct range
er = complex(max(1, real(er)), -abs(imag(er)));
ur = complex(max(1, real(ur)), -abs(imag(ur)));
thk = abs(thk);

%% Assemble Frequency Info (k_A1, k_A2, k_b1, k_b2)
k0(1, 1, 1, :) = 2*pi .* f(:) ./ O.c;
[k_A1, k_A2, k_b1, k_b2] = O.constructFrequencyMultipliers(k0);

%% Get A2, b2, and etaR1
A2 = O.A2;
b2 = O.b2;

etaR1 = sqrt(ur(1) ./ er(1));

%% Compute Required Tolerances
b = reshape(etaR1.*b2(1, 1).*k_b2(1, :, :), [], 1);
a = reshape(k_A1(1, 1, :), [], 1);

% AbsTolY = abs(a) ./ (1 + 2*abs(b).^2);

%% Check Error in Lossy Structure
[specE, specH] = O.multilayerSpectrumRect(O.lossy_tau, k0, er, ur, thk);

A1b1DomMode = sum(specE .* O.lossy_A1b1_E(:, 1, 1) ...
    + specH .* O.lossy_A1b1_H(:, 1, 1), 1);
errA1b1DomMode = sum(specE .* O.lossy_errA1b1_E(:, 1, 1) ...
    + specH .* O.lossy_errA1b1_H(:, 1, 1), 1);

% AbsTolY = AbsTol ./ abs(2 * a.*b ./ (a.*A1b1DomMode(:) + b).^2);

isLossyFrequency = (abs(errA1b1DomMode) ./ abs(A1b1DomMode)) < AbsTol;
% isLossyFrequency = abs(errA1b1DomMode(:)) < 0*AbsTolY(:);

%% Compute A1 and b1 (A1b1) at Lossy Frequencies
A1b1 = zeros(1, size(A2, 1), size(A2, 1) + 1, length(k0));

A1b1(:, :, :, isLossyFrequency) = ...
    sum(specE(:, 1, 1, isLossyFrequency) .* O.lossy_A1b1_E ...
    + specH(:, 1, 1, isLossyFrequency) .* O.lossy_A1b1_H, 1);

if O.verbose
    numLossyFreq = sum(isLossyFrequency);
    fprintf("Fixed-point integral method used for (%d/%d) frequencies.\n", ...
        numLossyFreq, length(f));
    if numLossyFreq < length(f)
        fprintf("Using adaptive integral method for the remaining (%d) frequencies.\n", ...
            length(f) - numLossyFreq);
    end
end

%% Compute A1 and b1 (A1b1) at Low Loss Frequencies
for ff = 1:length(f)
    if ~isLossyFrequency(ff)
        if O.verbose
            fprintf("f = %.4g GHz - ", f(ff));
        end
        A1b1(:, :, :, ff) = O.integralVectorized(...
            @(tauP) O.integralFunction(tauP, k0(ff), er, ur, thk), ...
            0, 1, AbsTol, 1e-8, 1, O.verbose);
%         A1b1(:, :, :, ff) = O.integralVectorized(...
%             @(tauP) O.integralFunction(tauP, k0(ff), er, ur, thk), ...
%             0, 1, AbsTol, AbsTolY(ff), 1, O.verbose);
    end
end

%% Extract A1 and b1 from A1b1
A1 = permute(A1b1(1, :, 1:end - 1, :), [2, 3, 4, 1]);
b1 = permute(A1b1(1, :, end, :), [2, 3, 4, 1]);

%% Calculate Reflection Coefficient at each Frequency
gam = zeros(length(f), 1);
for ff = 1:length(f)
    x = (A1(:, :, ff).*k_A1(:, :, ff) + etaR1.*A2.*k_A2(:, :, ff)) ...
        \ (b1(:, :, ff).*k_b1(:, :, ff) + etaR1.*b2.*k_b2(:, :, ff));
    gam(ff) = x(1);
end

end

