clc;
clear;
close all;

%% Inputs
k(1, :) = linspace(0, 10, 101);
a = 1;
order = 1;
kc = 1.51;

%% Define Integrand
besselj_deriv = @(order, x) 0.5 * (besselj(order - 1, x) - besselj(order + 1, x));

for ii = 1:numel(k)
    valNumeric(ii) = integral(...
        @(rho) rho .* besselj(order, kc .* rho) ...
        .* besselj(order, k(ii) .* rho), 0, a);
end

for ii = 1:numel(k)
    valNumericP(ii) = integral(...
        @(rho) rho .* besselj_deriv(order, kc .* rho) ...
        .* besselj(order, k(ii) .* rho), 0, a);
end

%% Compare
valTheory = a .* (k .* besselj(order, a.*kc) .* besselj(order - 1, a.*k) ...
    - kc .* besselj(order - 1, a.*kc) .* besselj(order, a.*k)) ...
    ./ (kc.^2 - k.^2);

valTheoryP = (1 - k) ./ (1 - k.*kc) ...
    .* (besselj(order, a.*kc) .* besselj(order, a.*k) ...
    - valTheory);

%% Plot
figure;
plot(k, valNumeric, "", LineWidth=1.5);
hold on;
plot(k, valTheory, "", LineWidth=1.5);
grid on;
ylim([min(valNumeric), max(valNumeric)]);

figure;
plot(k, valNumericP, "", LineWidth=1.5);
hold on;
plot(k, valTheoryP, "", LineWidth=1.5);
grid on;
ylim([min(valNumericP), max(valNumericP)]);

% figure;
% plot(k, valNumericP ./ valNumeric, "", LineWidth=1.5);
% hold on;
% % plot(k, valTheory, "", LineWidth=1.5);
% grid on;
% % ylim([min(valNumericP), max(valNumericP)]);



