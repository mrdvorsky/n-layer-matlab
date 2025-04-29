clc;
clear;
close all;

%% Inputs
dx = 0.125;
dy = dx;
k = 2*pi;
z0 = -5;

% Array Coordinates
xa(:, 1) = -20:dx:20;
ya(1, :) = -20:dx:20;

% Plot Coordinates
xp(:, 1) = -2:0.05:2;
yp(1, :) = -2:0.05:2;

%% Define Pattern Functions
coord_theta = @(x0, y0, z0) acos(-y0 ./ hypotn(x0, y0, z0));
coord_phi = @(x0, y0, z0) asin(x0 ./ hypot(x0, z0));

% Pattern for dipole with y-polarization.
pat_mag = @(x0, y0, z0) cos(coord_phi(x0, y0, z0)) .* exp(-(hypot(x0, y0)./abs(z0)).^2)...
    .* sin(coord_theta(x0, y0, z0)) ./ hypotn(x0, y0, z0);
pol_Ex = @(x0, y0, z0) cos(coord_theta(x0, y0, z0)) .* sin(coord_phi(x0, y0, z0));
pol_Ey = @(x0, y0, z0) sin(coord_theta(x0, y0, z0));
pol_Ez = @(x0, y0, z0) -cos(coord_theta(x0, y0, z0)) .* cos(coord_phi(x0, y0, z0));

pol_Hx = @(x0, y0, z0) cos(coord_phi(x0, y0, z0));
pol_Hy = @(x0, y0, z0) sin(coord_phi(x0, y0, z0)) .* cos(coord_theta(x0, y0, z0));
pol_Hz = @(x0, y0, z0) sin(coord_phi(x0, y0, z0)) .* sin(coord_theta(x0, y0, z0));

pat_Ex = @(x0, y0, z0) pat_mag(x0, y0, z0) .* pol_Ex(x0, y0, z0);
pat_Ey = @(x0, y0, z0) pat_mag(x0, y0, z0) .* pol_Ey(x0, y0, z0);
pat_Ez = @(x0, y0, z0) pat_mag(x0, y0, z0) .* pol_Ez(x0, y0, z0);

pat_Hx = @(x0, y0, z0) pat_mag(x0, y0, z0) .* pol_Hx(x0, y0, z0);
pat_Hy = @(x0, y0, z0) pat_mag(x0, y0, z0) .* pol_Hy(x0, y0, z0);
pat_Hz = @(x0, y0, z0) pat_mag(x0, y0, z0) .* pol_Hz(x0, y0, z0);

%% Calculate Array Coefficients
R = hypotn(xa, ya, z0);
S_meas = exp(1j .* k .* R) ./ R;

%% Calculate Focused Beam Pattern
% xa = gpuArray(xa);
% ya = gpuArray(ya);
% S_meas = gpuArray(S_meas);
% for xx = flip(1:numel(xp))
%     xx
%     for yy = flip(1:numel(yp))
%         R = hypotn(xp(xx) - xa, yp(yy) - ya, z0);
%         Ex(xx, yy) = gather(sum(S_meas .* exp(-1j .* k .* R) ...
%             .* pat_Ex(xp(xx) - xa, yp(yy) - ya, z0), ...
%             [1, 2]));
%         Ey(xx, yy) = gather(sum(S_meas .* exp(-1j .* k .* R) ...
%             .* pat_Ey(xp(xx) - xa, yp(yy) - ya, z0), ...
%             [1, 2]));
%         Ez(xx, yy) = gather(sum(S_meas .* exp(-1j .* k .* R) ...
%             .* pat_Ez(xp(xx) - xa, yp(yy) - ya, z0), ...
%             [1, 2]));
%     end
% end

%% Plotting
m = 25;
figure;
showImage(xa, ya, m .* pat_Hx(xa, ya, z0));

figure;
showImage(xa, ya, m .* pat_Hy(xa, ya, z0));

figure;
showImage(xa, ya, m .* pat_Hz(xa, ya, z0));


%% Plotting Beam
m = 1 ./ max(abs(Ey), [], "all");
figure;
showImage(xp, yp, m .* Ex, DisplayFormat="Magnitude");
clim([0, 1]);

figure;
showImage(xp, yp, m .* Ey, DisplayFormat="Magnitude");
clim([0, 1]);

figure;
showImage(xp, yp, m .* Ez, DisplayFormat="Magnitude");
clim([0, 1]);


