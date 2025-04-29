clc;
clear;
close all;

%% Inputs
krSamp = 10 * linspace(-1, 1, 1000);

fun1 = @(kr) sinc(kr).^2;

krR(:, 1) = 20 * linspace(-1, 1, 500);
krI(1, :) = 5 * linspace(-1, 1, 500);

%% Function
F = fun1(krSamp);
x(1, 1, :) = fftCoordinates(krSamp, ApplyFftShift=true);
f(1, 1, :) = fftshift(ifft(ifftshift(F)));

fun2 = @(kr) pi * sum(exp(-1j .* x .* kr) .* f, 3) .* abs(x(2) - x(1));


%% Plotting
figure;
plots(x, real(f), "", LineWidth=1.5);
grid on;

figure;
showImage(krR, krI, fun1(krR + 1j*krI), DisplayFormat="dB", Normalize=false);
axis normal;
clim([-150, 50]);

figure;
showImage(krR, krI, fun2(krR + 1j*krI), DisplayFormat="dB", Normalize=false);
axis normal;
clim([-150, 50]);

% figure;
% plot(krR, fun1(krR), "", LineWidth=1.5);
% hold on;
% plot(krR, fun2(krR), "", LineWidth=1.5);
% grid on;



