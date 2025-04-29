clc;
clear;
close all;

%% Inputs
nx = 100000;
xk = 2000;
L = 3000;
s = 2000;

%% Create Inputs
x(:, 1) = linspace(0, 1, nx).^2 + (0.5./nx)*rand(1, nx);
k(:, 1) = L./xk * (1 - linspace(1, 1./xk, xk)) ./ linspace(1, 1./xk, xk);

y = sqrt(s ./ pi) * exp(-s*(x - 0.5).^2);

%% Compute Spectrum
y_spec_theory = exp(-pi.^2 .* k.^2 ./ s) .* exp(-1j .* 2*pi*k .* 0.5);

tic;
y_spec1 = nufft(y, x, k) ./ nx;
toc;


y_spec2 = nufft_type3(y, x, k);


%% Plotting
% figure;
% plot(x, y, "", LineWidth=1.5);
% grid on;
% xlabel("x");
% ylabel("y");
% xlim([0, 1]);

% figure;
% loglog(k, abs(y_spec1), "o", LineWidth=1.5);
% hold on;
% loglog(k, abs(y_spec2), "x", LineWidth=1.5);
% loglog(k, abs(y_spec_theory), "", LineWidth=1.5);
% grid on;
% xlabel("k");
% ylabel("Y");
% ylim([1e-10, 1]);

figure;
loglog(k, abs(y_spec1 - y_spec_theory), "o", LineWidth=1.5);
hold on;
loglog(k, abs(y_spec2 - y_spec_theory), "x", LineWidth=1.5);
loglog(k, abs(y_spec2 - y_spec1), ".", LineWidth=1.5);
% loglog(k, abs(y_spec_theory), "", LineWidth=1.5);
grid on;
xlabel("k");
ylabel("error in Y");
ylim([1e-10, 1]);


