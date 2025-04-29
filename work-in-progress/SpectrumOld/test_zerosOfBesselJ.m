clc;
clear;
close all;

%% Inputs
% a(1, :) = [0.1, 0.25, 0.5, 0.75, 0.9];
a(1, :) = [0.1];
v = 1000;

r(:, 1) = linspace(0, 1100, 100001);

%% Calculate
% y = besselj(v, r);
y = besselj_derivative(v, r);

%% Plotting
figure;
for aa = 1:numel(a)
    [~, peakInds] = findpeaks(-y(:, aa).^2);
    bzeros = r(peakInds);
    plot(r(peakInds), "o", LineWidth=1.5);
    hold on;

    pFit = polyfit(1:numel(bzeros), bzeros, 1);
    plot(0:numel(bzeros), polyval(pFit, 0:numel(bzeros)), ":", LineWidth=1.5);


    for ii = 1:numel(bzeros)
        tic;
%         bzeros2(ii) = besselj_zeros(v, ii);
        bzeros2(ii) = besselj_derivative_zeros(v, ii);
        toc;
    end

    plot(bzeros2, "x", LineWidth=1.5);
end

bzeros(1)
bzeros2(1)

figure;
plot(r, db(y.^2), "", LineWidth=1.5);
grid on;
legend(compose("a = %g", a));


%% Functions
function [x] = besselj_zeros(v, n)
    fun = @(y) besselj(v, y);
    x = fzero(fun, v + 2.41*(v < 10));
    for ii = 2:n
        x_guess = x + pi*[0.9, 1.1];
        while sum(sign(fun(x_guess))) ~= 0
            x_guess(2) = x + (x_guess(2) - x)*1.1;
        end
        x = fzero(fun, x_guess);
    end
end

function [y] = besselj_derivative(v, x)
    y = 0.5 * (besselj(v - 1, x) - besselj(v + 1, x));
end

function [x] = besselj_derivative_zeros(v, n)
    fun = @(y) besselj_derivative(v, y);
    x = fzero(fun, v + 2.4*(v < 10));
    for ii = 2:n
        x_guess = x + pi*[0.9, 1.1];
        while sum(sign(fun(x_guess))) ~= 0
            x_guess(2) = x + (x_guess(2) - x)*1.1;
        end
        x = fzero(fun, x_guess);
    end
end


