clc;
clear;
close all;

%% Inputs
f = @(xs, ys) 1e1 * (1 + xs.^2 + 1.1*ys.^2)*0 + 1;

kPhi = 30;
kRho(:, 1) = linspace(0, 100, 1001);

Tx(:, 1) = 1.0 * [0, 1, 0] + 0;
Ty(:, 1) = 1.0 * [0, 0, 1] + 0;

x(:, 1) = linspace(min(Tx) - 0.1, max(Tx) + 0.1, 2001);
y(1, :) = linspace(min(Ty) - 0.1, max(Ty) + 0.1, 2001);

%% Sample Function
f_cropped = @(xs, ys) f(xs, ys) .* inpolygon(xs + 0*ys, 0*xs + ys, Tx, Ty);
f_sampled = f_cropped(x, y);

%% Perform Integrals
for ii = 1:numel(kRho)
    kx = kRho(ii) .* cosd(kPhi);
    ky = kRho(ii) .* sind(kPhi);

    g = @(xs, ys) f(xs, ys) .* exp(1j .* (kx.*xs + ky.*ys));
    g_cropped = @(xs, ys) f_cropped(xs, ys) .* exp(1j .* (kx.*xs + ky.*ys));

    int2(ii) = TriIntegral(g, Tx, Ty);
    disp(ii)
end

%% Plotting
% figure;
% surf(x, y, f_sampled.', LineStyle="none");
% hold on;
% plot3(Tx, Ty, f(Tx, Ty), "kx", LineWidth=1.5);
% plot3(Tx_mid, Ty_mid, f(Tx_mid, Ty_mid), "ko", LineWidth=1.5);

figure;
plot(kRho, abs(int2), "", LineWidth=1.5);
grid on;

figure;
plot(kRho, rad2deg(angle(int2)), "", LineWidth=1.5);
grid on;







