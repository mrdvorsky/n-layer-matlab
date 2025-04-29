clc;
clear;
close all;

%% Inputs
p = [0, 1];
x1 = [-0.5];
x2 = [1];

% p = [2, -3, 1, -6];
% x1 = [3, 2];
% x2 = [2, 0];

kr = fejer2_halfOpen(5000, 1);

%% Quadrature
[xs, ws] = fejer2(2000, min([x1, x2]), max([x1, x2]));

spec1 = nufft(polyval(p, xs) .* ws, xs, kr ./ (2*pi), 1);

%% Analytical
[x0, s0, c0] = polyFourierTransform2(x1, x2, ...
    [polyval(p, linspace(x1(1), x2(1), numel(p)))]);
% [x0, s0, c0] = polyFourierTransform2(x1, x2, ...
%     [polyval(p, linspace(x1(1), x2(1), numel(p))); ...
%     polyval(p, linspace(x1(2), x2(2), numel(p)))]);

spec2 = 0;
for nn = 1:numel(p)
    spec2 = spec2 + sum(c0(:, nn).' .* exp(-1j .* x0(:).' .* kr)...
        .* besselj_spherical(nn - 1, s0(:).' .* kr), 2);
end


%% Plot
figure;
plot(1:numel(kr), real(spec1), "", LineWidth=1.5);
hold on;
plot(1:numel(kr), imag(spec1), "", LineWidth=1.5);
grid on;

figure;
plot(1:numel(kr), real(spec2), "", LineWidth=1.5);
hold on;
plot(1:numel(kr), imag(spec2), "", LineWidth=1.5);
grid on;









