clc;
clear;
% close all;

%% Inputs
kx(:, 1) = linspace(-0.5, 5, 1600);
ky(1, :) = linspace(-1.5, 1.5, 1200);

k0 = 1;
er = {4 - 0.001j};
ur = {1};
thk = {0.001};

% numPoles = 5;

%% Calculate
% [GamH, GamE] = nLayer.computeGamma0(kx + 1j*ky, k0, er, ur, thk);

gamFun = @(x) nLayer.computeGammaA(x, k0, er, ur, thk);
GamH = gamFun(kx + 1j*ky);

%% Pole Locations
% m = 1:numPoles;
% kr_poles = sqrt(er{1}.*k0.^2 - (m * pi ./ thk{1}).^2);
% 
% for ii = 1:numel(m)
%     path = kr_poles(ii) + [1 + 1j, -1 + 1j, -1 - 1j, 1 - 1j];
%     res(ii) = integral(gamFun, path(1), path(1), Waypoints=path);
% end
% res



%% Plot
figure(1);
hold off;
showImage(kx, ky, GamH.^30, DisplayFormat="MagPhase");
grid on;
% clim(1*[0, 1]);
hold on;
% plot(kr_poles, "xw");





