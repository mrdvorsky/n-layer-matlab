clc;
clear;
close all;

%% Inputs
w = @(x) 100*(besselj(2, 0.1 * x) + 1j*besselj(2, 0.1 * x)) .* (1 + besselj(0, 100*x)) .* exp(-2.1 * x) .* x.^2;
f = @(x) polyval([1, 1, -1, 1], x) ./ (1 + x.^3);

a = 0.0001;
b = 1.1;

intOrder = 10;

L = 6;

%% Integrate
I1 = integral(@(x) w(x) .* f(x), 0, inf);

[x2, x2_weights] = clenshawCurtisHalfOpen(intOrder, L);
I2 = sum(w(x2) .* f(x2) .* x2_weights);

[x3, x3_weights] = clenshawCurtisHalfOpen(intOrder, L, WeightingFunction=w);
I3 = sum(f(x3) .* x3_weights);

%% Error
err2_db = db(I1 - I2)
err3_db = db(I1 - I3)







