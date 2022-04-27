function [k_A1, k_A2, k_b1, k_b2] = constructFrequencyMultipliers(O, f)
%CONSTRUCTFREQUENCYMULTIPLIERS Computes k_A1, k_A2, k_b1, k_b2.
% Computes the matrices k_A1, k_A2, k_b1, k_b2, used to solve the 
% rectangular waveguide equation for multilayer structures.
%
% Inputs:
%   f - vector of frequencies to consider (default unit is GHz).
% Outputs:
%   k_A1, k_A2, k_b1, k_b2 - Matrices used to compute S11. See usage below.
%
% The outputs of this function can be used along with the outputs of the
% constructMatrixEquation(...) function to calculate S11 for a
% rectangular waveguide. See example usage below, where "f" is scalar.
%
% Example Usage:
%   [k_A1, k_A2, k_b1, k_b2] = O.constructFrequencyMultipliers(f);
%   [A1, A2] = O.constructMatrixEquation(nLayerInt);
%   etaR1 = sqrt(ur(1) ./ er(1));
%   x = (A1.*k_A1 + etaR1.*A2.*k_A2) ...
%        \ (-A1(:, 1).*k_b1 + etaR1.*A2(:, 1).*k_b2);
%   S11 = x(1);
%
% Although the example above shows usage with a scalar value for "f", the
% input "f" can be a vector. In this case, the size of the 3rd dimension
% of each output matrix will be equal to numel(f).
%
% Author: Matt Dvorsky

arguments
    O;
    f(1, 1, :) double;
end

%% Mode Coefficients
k0 = 2*pi .* f ./ O.speedOfLight;

am(1, :) = O.modesTE(:, 1) * pi ./ O.a;
bn(1, :) = O.modesTE(:, 2) * pi ./ O.b;
ap(:, 1) = O.modesTE(:, 1) * pi ./ O.a;
bq(:, 1) = O.modesTE(:, 2) * pi ./ O.b;

kmn = conj(sqrt(k0.^2 - am.^2 - bn.^2));
kpq = conj(sqrt(k0.^2 - ap.^2 - bq.^2));

indTE = find(O.modesTE(:, 1) > 0);
indTM = find(O.modesTE(:, 2) > 0);

%% Assemble Frequency Multipliers
k_A1 = [k0 + 0*kmn(1, indTE, :), 0*k0 + kmn(1, indTM, :)];
k_A2 = [0*k0 + kmn(1, indTE, :), k0 + 0*kmn(1, indTM, :)];

k_b1 = [k0 + 0*kpq(indTE, 1, :); k0 + 0*kpq(indTM, 1, :)];
k_b2 = [0*k0 + kpq(indTE, 1, :); 0*k0 + kpq(indTM, 1, :)];

end

