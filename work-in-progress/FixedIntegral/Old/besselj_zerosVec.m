function [jvm] = besselj_zerosVec(v, n)
%BESSELJ_ZEROSVEC Zeros of bessel function.
%
% Author: Matt Dvorsky

arguments
    v(1, 1) {mustBeNonnegative};
    n(1, 1) {mustBePositive, mustBeInteger};
end

%% Get Initial Guess
jvm(:, 1) = ((1:n) + 0.5*v - 0.25) * pi;

%% Refine Initial Guess
for ii = 1:numel(jvm)
    jvm(ii) = fzero(@(x) besselj(v, x), jvm(ii) + pi*[-0.5, 0.5]);
end

end

