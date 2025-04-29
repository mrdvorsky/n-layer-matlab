clc;
clear;
close all;

%% Inputs
N = 20;

Nm = 75;

%% Convergence
Nrho = 400;
Nphi = 150;

L = 1;
Lc = 1;

wgA = 1;
wgB = 0.5;

%% Helper
[kr(:, 1), kr_weights(:, 1)] = fejer2_halfOpen(Nrho, L);
[krc, momentKernel] = getContourIntegrals(N, L, Lc);


%% Weight Functions
% for kk = 1:Nrho
%     kk
%     w = @(x) momentKernel(x, Nm) ./ 1e12;
%     M_k(kk) = 1e12 * integral(@(x) (2*L) * w(L * cot(0.5*x).^2) ...
%         .* sin(kk*x), 0, pi);
% end

w = @(x) 4 * momentKernel(x, Nm) ./ 1e12 ./ (1 + x).^2;
[~, kr_weights2(:, 1)] = fejer2_halfOpen(Nrho, L, WeightingFunction=w);
kr_weights2 = kr_weights2 .* (1 + kr).^2;

%% Plotting
figure;
plot(1:numel(kr), real(momentKernel(kr, Nm)), "", LineWidth=1.5);
hold on;
plot(1:numel(kr), imag(momentKernel(kr, Nm)), "", LineWidth=1.5);
grid on;

figure;
plot(1:numel(kr), real(kr_weights2), "", LineWidth=1.5);
hold on;
plot(1:numel(kr), imag(kr_weights2), "", LineWidth=1.5);
grid on;


%%
krVal = 0.4746;

for ii = 1:Nm
    kVal(ii, 1) = momentKernel(krVal, ii);
end

kValSpec = fft([0; kVal; 0; -flip(kVal)]);

figure;
plot(1:numel(kVal), real(kVal), "", LineWidth=1.5);
hold on;
plot(1:numel(kVal), imag(kVal), "", LineWidth=1.5);
grid on;

figure;
plot(1:numel(kValSpec), real(kValSpec), "", LineWidth=1.5);
hold on;
plot(1:numel(kValSpec), imag(kValSpec), "", LineWidth=1.5);
grid on;










