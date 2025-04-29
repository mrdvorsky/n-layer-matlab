clc;
clear;
close all;

%% Inputs
kx(:, 1) = linspace(-0.5, 5, 600);
ky(1, :) = linspace(-1.5, 1.5, 200);

% krc1(:, 1) = fejer2(200, 1j, 10 + 1j);
% krc2(:, 1) = fejer2(200, 0.8j, 10 + 0.8j);

krc1(:, 1) = fejer2_halfOpen(200, 0.3) + 2j;
krc2(:, 1) = fejer2_halfOpen(200, 0.3) + 1.8j;

krc = ([krc1; krc2]);

f_res = 37.474;
k0 = 2*pi*f_res / 299.792468;

Lcw = 20;
Lch = 1;

Nm = 5;

er = {1 - 0.1j};
ur = {1};
thk = {22};

%% Calculate 
[GamH, GamE] = nLayer.computeGamma0(kx + 1j*ky, k0, er, ur, thk);
[GamH_c, GamE_c] = nLayer.computeGamma0(krc, k0, er, ur, thk);

%% Plot
figure(1);
subplot(2, 2, 1);
ImgH = showImage(kx, ky, GamH, DisplayFormat="Magnitude");
grid on;
clim(10*[0, 1]);
% colormap colormapplusminus;

subplot(2, 2, 2);
ImgE = showImage(kx, ky, GamE, DisplayFormat="Magnitude");
grid on;
clim(10*[0, 1]);
% colormap colormapplusminus;

%% Plot Integral Paths
subplot(2, 2, 3);
contourH = plot([real(GamH_c), imag(GamH_c)], "", LineWidth=1.5);
grid on;

subplot(2, 2, 4);
contourE = plot([real(GamE_c), imag(GamE_c)], "", LineWidth=1.5);
grid on;

%% Interactive
subplot(2, 2, 1);
hold on;
interactivePlot([0.4, 1], [-0.2, -0.1], ...
    {@updatePlot, ImgH, ImgE, contourH, contourE, kx, ky, krc, k0}, ...
    DragClampFun=@(x, y) [max(0, x), min(0, y)], ...
    MarkerSize=20, MarkerColor="w");



%% Helper
function [x, y] = updatePlot(x, y, ind, ImgH, ImgE, contourH, contourE, kx, ky, krc, k0)
    kr = x + 1j*y;
    phase_val = 2*2*pi + pi;
    
    erp = (real(kr(2)).^2 - imag(kr(2)).^2) ./ k0.^2;
    erpp = -imag((kr(2) ./ k0).^2);

    if (ind == 1)
        erpp = -imag(kr(1).^2 ./ k0.^2);
    end
    
    kz = sqrt(k0.^2 .* (erp - 1j*erpp) - kr(1).^2);
    thk = { phase_val ./ real(2*kz) };
    er = { erp - 1j*erpp };

    kr(2) = sqrt(er{1}) .* k0;
    kr(1) = -1j*sqrt((0.5*phase_val./thk{1}).^2 - er{1}.*k0.^2);

    x = real(kr);
    y = imag(kr);


    % Update Images
    [GamH, GamE] = nLayer.computeGamma0(kx + 1j*ky, k0, er, {1}, thk);
    % [GamH, GamE] = nLayer.computeGamma0(kx + 1j*ky, k0, er, {1}, thk);
    ImgH.CData = (real(GamE)).';
    ImgE.CData = (imag(GamH)).';
    title(ImgH.Parent, sprintf("t = %g mm, er = %g - j%g", ...
        thk{1}, real(er{1}), imag(er{1})));

    % Update Plots
    [GamH_c, GamE_c] = nLayer.computeGamma0(krc, k0, er, {1}, thk);
    % [GamH_c, GamE_c] = nLayer.computeGamma0(krc, k0, er, {1}, thk);
    contourH(1).YData = real(GamH_c ./ krc);
    contourH(2).YData = imag(GamH_c ./ krc);
    contourE(1).YData = real(GamE_c .* krc);
    contourE(2).YData = imag(GamE_c .* krc);
end








