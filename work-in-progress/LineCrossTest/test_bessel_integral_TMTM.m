clc;
clear;
close all;

%% Inputs
kc = 1;
kr(1, :) = [linspace(0.0001, 0.1, 50), linspace(0.1, 10, 400)];

L = 2;
M = 10;

%% Nodes
[x, w] = fejer2_halfOpen(100000, 1);

fun = @(x, kr) besselj(0, kr.*x) .* x.^2 ./ (x.^2 - kc.^2).^2;

fun_m = @(x, kr, m) fun(x, kr) .* sin(2*m * acot(sqrt((x - 2)./L))) ...
                               ./ sin(2   * acot(sqrt((x - 2)./L)));

for ii = 1:numel(kr)
    val(ii) = integral(@(x) fun_m(x, kr(ii), M), 2, 50000);
end

% for ii = 1:numel(kr)
%     val(ii) = integral(@(x) fun_m(x, kr(ii), M), 1, 1 + 0.1j) ...
%         + integral(@(x) fun_m(x, kr(ii), M), 1 + 0.1j, 2) ...
%         + integral(@(x) fun_m(x, kr(ii), M), 2, 5000);
% end

%% Plotting
figure;
plot(kr, real(val));
hold on;
plot(kr, imag(val));

% figure;
% plot(kr, abs(val) ./ (1 + kr));
% hold on;
% plot(kr, imag(val));
% ylim([-1, 5])


