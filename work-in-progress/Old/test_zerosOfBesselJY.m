clc;
clear;
close all;

%% Inputs
a(1, :) = [0.5];
v = 0;

r(:, 1) = linspace(0, 30, 100001) ./ (1 - a);

%% Calculate
y = besseljy_det(v, r, a);
% y = besseljy_derivative_det(v, r, a);

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
        bzeros2(ii) = besseljy_det_zeros(v, ii, a);
%         bzeros2(ii) = besseljy_det_derivative_zeros(v, ii, a);
        toc;
    end

    plot(bzeros2, "x", LineWidth=1.5);
end

bzeros(1)
bzeros2(1)

figure;
plot(r, db(y), "", LineWidth=1.5);
grid on;
legend(compose("a = %g", a));


%% Functions
function [x] = besseljy_det_zeros(v, n, a)
    fun = @(y) besseljy_det(v, abs(y), a);
    x = abs(fzero(fun, (v + (v == 0))));
    for ii = 2:n
        x_guess = x + (0.5*pi ./ (1 - a))*[0.9, 1.1];
        while sum(sign(fun(x_guess))) ~= 0
            x_guess(2) = x + (x_guess(2) - x)*1.1;
        end
        x = fzero(fun, x_guess);
    end
end

function [x] = besseljy_det_derivative_zeros(v, n, a)
    fun = @(y) besseljy_derivative_det(v, abs(y), a);
    x = abs(fzero(fun, (v + (v == 0))));
    for ii = 2:n
        x_guess = x + (0.5*pi ./ (1 - a))*[0.9, 1.1];
        while sum(sign(fun(x_guess))) ~= 0
            x_guess(2) = x + (x_guess(2) - x)*1.1;
        end
        x = fzero(fun, x_guess);
    end
end

function [y] = besseljy_det(v, x, a)
    y = besselj(v, x).*bessely(v, a.*x) ...
        - besselj(v, a.*x).*bessely(v, x);
end

function [y] = besseljy_derivative_det(v, x, a)
    y = 0.5 * (besseljy_det(v - 1, x, a) - besseljy_det(v + 1, x, a));
end



