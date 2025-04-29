clc;
clear;
close all;

%% Inputs
w = @(r) exp(2j * r) ./ (1 + r).^2;

n = 1;
L = 3;

%% Compute Integrals
I1 = integral(@(x) w(L * cot(0.5*x).^2) ...
            .* sin(x) .* cos(n*x), ...
            0, pi);

test2 = @(x) 2 * w(x) ./ (1 + x).^2 .* cos(2 * n .* acot(sqrt(x)));
I2 = integral(test2, 0, inf);

[x3, x3_weights] = clenshawCurtisHalfOpen(200, 35);
I3 = sum(test2(x3) .* x3_weights);

%% Error
err2_db = db(I1 - I2)
err3_db = db(I1 - I3)




