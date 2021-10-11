function [k_A1, k_A2, k_b1, k_b2] = constructFrequencyMultipliers(O, k0)
%CONSTRUCTFREQUENCYMULTIPLIERS Summary of this function goes here
%   Detailed explanation goes here

arguments
    O
    k0(1, 1, :) {mustBeNumeric}
end

%% Mode Coefficients
am(1, :) = O.modes(:, 1) * pi ./ O.a;
bn(1, :) = O.modes(:, 2) * pi ./ O.b;
ap(:, 1) = O.modes(:, 1) * pi ./ O.a;
bq(:, 1) = O.modes(:, 2) * pi ./ O.b;

kmn = conj(sqrt(k0.^2 - am.^2 - bn.^2));
kpq = conj(sqrt(k0.^2 - ap.^2 - bq.^2));

indTE = find(O.modes(:, 1) > 0);
indTM = find(O.modes(:, 2) > 0);

%% Assemble Frequency Multipliers
k_A1 = [k0 + 0*kmn(1, indTE, :), 0*k0 + kmn(1, indTM, :)];
k_A2 = [0*k0 + kmn(1, indTE, :), k0 + 0*kmn(1, indTM, :)];

k_b1 = [k0 + 0*kpq(indTE, 1, :); k0 + 0*kpq(indTM, 1, :)];
k_b2 = [0*k0 + kpq(indTE, 1, :); 0*k0 + kpq(indTM, 1, :)];

end

