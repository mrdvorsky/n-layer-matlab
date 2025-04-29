clc;
clear;
close all;

%% Inputs
N = 3;
L = 1;
Lc = 1.0;

Nrho = 16000;
Nphi = 1*451;

%% Structure
er = 2 - 0.0001j;
ur = 1;
thk = 2.5;
f = 30;
k0 = f * (2*pi/299.792458);

%% nLayer
wgA = 1;
wgB = 0.5;

%% Get Spectrum Functions
[specX1, specY1] = nLayerOpenEnded.getSpectrumRectangular(...
    wgA, wgB, 3, 2, "TE");

%% Helper
% [phi1(1, 1, :)] = fejer2(Nphi, 0, pi/2);
[phi1(1, 1, :)] = linspace(0, pi, Nphi);
% phi1_weights = 4*phi1_weights;

modeFunE = @(kr, phi) kr .* sum(specY1(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* sin(phi).^2, 3);
modeFunH = @(kr, phi) kr .* sum(specY1(kr .* cos(phi), kr .* sin(phi)).^2 ...
    .* cos(phi).^2, 3);

%% Compute Weights and Nodes
[krc, momentKernel] = getContourIntegrals(N, L, Lc);
[kr(:, 1), kr_weights(:, 1)] = fejer2_halfOpen(Nrho, L);


for pp = 1:numel(phi1)
    pp
    MkE(pp, :) = sum((1 ./ sqrt(1 + kr.^2)) .* kr_weights .* modeFunE(kr, phi1(pp)) .* (momentKernel(kr, 1:N) - 0*(1:N)), 1);
    MkH(pp, :) = sum((sqrt(1 + kr.^2)) .* kr_weights .* modeFunH(kr, phi1(pp)) .* (momentKernel(kr, 1:N) - 0*(1:N)), 1);
end


%% Plotting
figure;
plots(rad2deg(phi1), imag(MkE), "", LineWidth=1.5);
grid on;

figure;
plots(rad2deg(phi1), imag(MkH), "", LineWidth=1.5);
grid on;

% [~, krc_weights_E(:, 1)] = fejer2_halfOpen(N, L, ...
%     WeightingMoments=MkE);
% 
% [~, krc_weights_H(:, 1)] = fejer2_halfOpen(N, L, ...
%     WeightingMoments=MkH);








