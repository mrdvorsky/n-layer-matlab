clc;
clear;
close all;

%% Inputs
N = 100;

k(:, 1) = 0:200;

%% Quad Rules
x0 = (0.5 + (0:N - 1)) ./ N;
w0 = 0*x0 + (1 ./ N);

[x1, w1] = fejer2(N, 0, 1);
[x2, w2] = gaussLegendre(N, 0, 1);

%% Integrate
for kk = 1:numel(k)
    kk

    % f = @(x) exp((2*pi) * 1j .* k(kk) .* x);

    f = @(x) exp(1j .* k(kk) .* cos(2*pi .* x));

    It(kk, 1) = integral(f, 0, 1);

    I0(kk, 1) = sum(f(x0) .* w0);
    I1(kk, 1) = sum(f(x1) .* w1);
    I2(kk, 1) = sum(f(x2) .* w2);
end


%% Plotting
figure;
plot(k, db(I0 - It), "", LineWidth=1.5);
hold on;
plot(k, db(I1 - It), "", LineWidth=1.5);
plot(k, db(I2 - It), "", LineWidth=1.5);
grid on;



