function [gam] = calculate(O, f, er, ur, thk, options)
%CALCULATE Calculate S11 for TE10 mode
% Inputs
%   f - vector of frequencies (GHz)
%   er - vector of complex relative permittivities for each layer
%   ur - vector of complex relative permeabilities for each layer
%   thk - vector of thicknesses for each layer (same length as er and ur)

arguments
    O;
    f(:, 1);
    er(1, :);
    ur(1, :);
    thk(1, :);
    options.AbsTol(1, 1) = O.convergenceAbsTol;
end

% Make sure er, ur, and thk are in the correct range
er = complex(max(1, real(er)), -abs(imag(er)));
ur = complex(max(1, real(ur)), -abs(imag(ur)));
thk = abs(thk);

%% Compute A1 and b1
% This is the computationally intensive part of this algorithm
[A1, b1] = O.computeA1b1(f, er, ur, thk, options.AbsTol);

%% Get A2, b2, and etaR1
% A2 and b2 are precomputed in the "recomputeInterpolants" function.
A2 = O.A2;
b2 = O.b2;

etaR1 = sqrt(ur(1) ./ er(1));

%% Assemble Frequency Info (k_A1, k_A2, k_b1, k_b2)
[k_A1, k_A2, k_b1, k_b2] = O.constructFrequencyMultipliers(f);

%% Calculate Reflection Coefficient at each Frequency
gam = zeros(length(f), 1);
for ff = 1:length(f)
    x = (A1(:, :, ff).*k_A1(:, :, ff) + etaR1.*A2.*k_A2(:, :, ff)) ...
        \ (b1(:, :, ff).*k_b1(:, :, ff) + etaR1.*b2.*k_b2(:, :, ff));
    gam(ff) = x(1);
end

end

