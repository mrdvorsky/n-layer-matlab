function [wFun0, wFun2] = getWeightingFunctions(krMax, rMin, rMax, Nr, Nkr)
%GETWEIGHTINGFUNCTIONS Summary of this function goes here
%   Detailed explanation goes here


r(:, 1) = linspace(0, rMax, Nr);

%% Integral
% [kr(1, :), kr_w(1, :)] = fejer2_halfOpen(Nkr, 1);
[kr(1, :), kr_w(1, :)] = fejer2(Nkr, 0, krMax);


val0 = pi * sum(kr ./ (1 + kr) .* besselj(0, kr .* r) .* kr_w, 2);
val2 = pi * sum(kr ./ (1 + kr) .* besselj(2, kr .* r) .* kr_w, 2);

% val0 = 1 ./ r;
% val2 = 1 ./ r;

wFun0 = griddedInterpolant(r, val0, "spline");
wFun2 = griddedInterpolant(r, val2, "spline");

end

