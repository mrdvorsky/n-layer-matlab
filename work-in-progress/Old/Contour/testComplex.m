clc;
clear;
close all;

%% Inputs
w = @(x) besselh(2, 100 * x) .* exp(-0.1 * x);
f = @(x) polyval([1, 1, -1, 1], x);

a = 0.001;
b = 1;

intOrder = 10;

%% Integrate
I1 = integral(@(x) w(x) .* f(x), a, b);

[x2, x2_weights] = clenshawCurtis(intOrder, a, b);
I2 = sum(w(x2) .* f(x2) .* x2_weights);

[x3, x3_weights] = clenshawCurtis(intOrder, a, b, WeightingFunction=w);
I3 = sum(f(x3) .* x3_weights);

%% Error
err2_db = db(I1 - I2)
err3_db = db(I1 - I3)







