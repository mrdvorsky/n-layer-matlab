clc;
clear;
% close all;

%% Inputs
krMax = 200;

Lc = 0.1;
ccL = 0.3;

ccOrder = 1;

r(:, 1) = linspace(-50, 50, 1000);

%% Helper
xToKrc = @(x) x + 1j .* Lc .* x ./ sqrt(Lc.^2 + x.^2);
xToKrc_weights = @(x) 1 + 1j .* Lc.^3 ./ (x.^2 + Lc.^2).^(1.5);

%% Integrate
for rr = 1:numel(r)
    rr
    w = @(x) (exp(-1j .* r(rr) .* xToKrc(x)) + exp(1j .* r(rr) .* xToKrc(x))) .* xToKrc_weights(x);

    val(rr) = integral(...
        @(x) w(x) ...
        .* (cos(2*ccOrder * acot(sqrt(x ./ ccL)))), ...
        0, krMax);
end

%% Plotting
figure;
plot(r, val, "", LineWidth=1.5);
grid on;



