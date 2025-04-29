clc;
clear;
close all;

%% Inputs
t1_x = [-1, 1, 0];
t1_y = [0, 0, 1];

t2_x = [-0.5, 0.5, 0];
t2_y = [0.5, 0.5, 1];

x(:, 1) = linspace(-2, 2, 1001);
y(1, :) = linspace(-2, 2, 1001);

%% Sample
E1 = 1.0 * inpolygon(x + 0*y, y + 0*x, t1_x, t1_y);
E2 = 1.0 * inpolygon(x + 0*y, y + 0*x, t2_x, t2_y);

E12 = fftshift(ifft2(fft2(E1) .* fft2(flip(flip(E2, 1), 2))));
E12 = abs(E12);
E12 = E12 ./ max(E12(:));

%% Plot
figure;
imgPlot = showImage(x, y, E12, DisplayFormat="Magnitude");
xlim([-2, 2]);
ylim([-2, 2]);
colormap parula;
hold on;
[~, contPlot] = contour(x, y, abs(E12).', "k");
interactivePlot([t1_x, t2_x], [t1_y, t2_y], ...
    {@updateImg, imgPlot, contPlot, x, y}, ...
    MarkerColor="k");




%% Helpers
function updateImg(xp, yp, imgPlot, contPlot, x, y)
    E1 = 1.0 * inpolygon(x + 0*y, y + 0*x, xp(1:3), yp(1:3));
    E2 = 1.0 * inpolygon(x + 0*y, y + 0*x, xp(4:6), yp(4:6));
    
    E12 = fftshift(ifft2(fft2(E1) .* fft2(flip(flip(E2, 1), 2))));
    E12 = abs(E12).';
    E12 = E12 ./ max(E12(:));
    
    imgPlot.CData = E12;
    contPlot.ZData = E12;
end






