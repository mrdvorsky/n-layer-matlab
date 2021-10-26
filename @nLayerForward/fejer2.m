function [nodes, weights, errorWeights] = fejer2(orderN, a, b)
%FEJER2 Generate Fejer Type II weights and nodes for closed interval integration.
% This function generates the weights and nodes required to compute a
% definite integral over a closed interval. The weights and nodes are
% defined using the Fejer Type II Quadrature rules. These rules are very
% similar to the Clenshaw-Curtis rules, but do not result in evalution of
% the funtion at the endpoints a and b.
%
% The function outputs "nodes" and "weights" can be used to approximate
% the definite integral of a function f(x)dx over the interval [a,b] by
% computing q = sum(weights .* f(nodes)). This should give approximately
% the same result as q = integral(f, a, b), with a higher value of
% orderN resulting in a better approximation. The error in q can be
% estimated using the output parameter errorWeights using the formula
% qErr = sum(errorWeights .* f(nodes)).
%
% If f(x) is a polynomial with degree less than orderN - 1, the result will
% be exact. Use of these quadrature rules will not result in evalation of
% the function at the interval endpoints a and b.
%
% Example Usage:
%   [nodes, weights] = fejer2(N, a, b);
%   [nodes, weights, errorWeights] = fejer2(N, a, b);
%   q = sum(fun(nodes) .* weights, 1);
%   qErr = sum(fun(nodes) .* errorWeights, 1);
%
% Inputs:
%   orderN - Scalar number of nodes to calculate.
%   a - Scalar integration lower bound. Must be finite.
%   b - Scalar integration upper bound. Must be finite.
% Outputs:
%   nodes - Column vector of coordinates at which to evaluate function.
%   weights - Column vector of weights to perform weighted sum.
%   errorWeights - Column vector of weights to estimate integration error.
%
% Author: Matt Dvorsky

arguments
   orderN(1, 1) {mustBeInteger, mustBeFinite} = 10;
    a(1, 1) {mustBeNumeric, mustBeFinite} = -1;
    b(1, 1) {mustBeNumeric, mustBeFinite} = 1;
end

%% Calculate Moments
N = orderN + 1;

if nargout == 3
    N = 2*round(0.5*N);
end

kBegin(:, 1) = 0:1:floor(0.5*N - 1);

momentsBegin = 2 ./ (1 - 4*kBegin.^2);
momentsMid = (N - 3) ./ (2*floor(0.5*N) - 1) - 1;

%% Calculate Weights
if (mod(N, 2) == 0) % If N is even
    weights = real(ifft([momentsBegin; momentsMid; flip(momentsBegin(2:end))]));
else
    weights = real(ifft([momentsBegin; momentsMid; momentsMid; flip(momentsBegin(2:end))]));
end

weights = weights(2:end, 1);

%% Calculate Nodes
nodes(:, 1) = -cos(pi * (1:N - 1) ./ (N));

%% Change Interval
weights = 0.5*(b - a) .* weights;
nodes = 0.5*(b - a) .* nodes + 0.5*(a + b);

%% Compute Error Estimate Weights
if nargout == 3
    [~, weightsReducedOrder] = ...
        fejer2(floor(0.5*N - 1), a, b);
    
    errorWeights = weights;
    errorWeights(2:2:end) = errorWeights(2:2:end) - weightsReducedOrder;
end

end

