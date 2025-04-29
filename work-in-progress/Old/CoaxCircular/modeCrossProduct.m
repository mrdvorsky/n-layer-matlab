function [A, coaxBeta, circBeta] = modeCrossProduct(coaxInnerR, coaxOuterR, circR, nCoax, nCirc)
%MODECROSSPRODUCT Summary of this function goes here.
% 
% Author: Matt Dvorsky

arguments
    coaxInnerR(1, 1);
    coaxOuterR(1, 1);
    circR(1, 1);
    nCoax(1, 1) {mustBePositive, mustBeInteger};
    nCirc(1, 1) {mustBePositive, mustBeInteger};
end

%% Calculate Mode Properties
for n = 1:nCoax
    [coaxBeta(n), coaxA(n), coaxB(n)] = ...
        coaxial_cutoffs_TM(coaxInnerR, coaxOuterR, 0, n);
end

for n = 1:nCirc
    [circBeta(n)] = ...
        circular_cutoffs_TM(circInnerR, 0, n);
end

%% Calculate Mode Fields




end


%% Helper Functions
function [y] = besseljy(a, b, v, x)
    y = a.*besselj(v, x) + b.*bessely(v, x);
end

function [y] = besseljy_derivative(a, b, v, x)
    y = 0.5 * (besseljy(a, b, v - 1, x) - besseljy(a, b, v + 1, x));
end

function [y] = besselj_derivative(v, x)
    y = 0.5 * (besselj(v - 1, x) - besselj(v + 1, x));
end

function [kc, a, b] = coaxial_cutoffs_TM(wgR1, wgR2, m, n)
    kc = besseljy_det_zeros(m, n, wgR1 ./ wgR2) ./ wgR2;
    ab = null([besseljy(1, 0, m, kc .* wgR2), ...
        besseljy(0, 1, m, kc .* wgR2)]);
    a = ab(1);
    b = ab(2);
end

function [kc] = circular_cutoffs_TM(wgR, m, n)
    kc = besselj_zeros(m, n) ./ wgR;
end

%% Bessel Function Zeros
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


