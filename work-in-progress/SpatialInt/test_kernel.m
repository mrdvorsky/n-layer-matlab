clc;
clear;
close all;

%% Inputs
Nm = 4;

x = linspace(0.0, 10, 1000);

Lch = 2;
Lcw = 10;

%% Helper
krToKrc = @(x) x + 1j * (Lch.*Lcw) * x ./ sqrt((Lch + x.^2).*(3*Lcw.^2 + x.^2));

%% Integrate
kern = @(kr) sin(2*(Nm) .* acot(sqrt(kr))) ...
          ./ sin(2 *       acot(sqrt(kr))) ...
          - Nm + 40 ./ (1 + kr) - 96 ./ (1 + kr.^2);

for xx = 1:numel(x)
    xx
    y(xx) = integral(@(k) besselj(0, x(xx) * k) .* k .* kern((k)), 0, 1000);
end

%% Plot
figure;
plot(x, real(y), "", LineWidth=1.5);
hold on;
plot(x, imag(y), ":", LineWidth=1.5);
grid on;





