function [nodes, weights] = clenshawCurtiss(orderN, a, b, options)
%CLENSHAWCURTIS Generate Clenshaw-Curtiss weights and nodes for closed interval integration.
% This function generates the weights and nodes required to compute a
% definite integral over a closed interval. The weights and nodes are
% defined using the Clenshaw-Curtis Quadrature rules.
%
% The function outputs "nodes" and "weights" can be used to approximate
% the definite integral of a function f(x)dx over the interval [a,b] by
% computing q = sum(weights .* f(nodes)). This should give approximately
% the same result as q = integral(f, a, b), with a higher value of
% orderN resulting in a better approximation.
%
% If f(x) is a polynomial with degree less than orderN, the result will be
% exact. Use of these quadrature rules will result in evalation of the
% function at the interval endpoints a and b.
%
% A weighting function "w" can be optionally supplied such that the
% integral 1 = sum(weights .* f(nodes)) corresponds to the integral of
% f(x)w(x)dx over the closed interval [a,b]. In this case, the value 1
% will be exact f(x) is a polynomial with degree less than orderN,
% regardless of w(x). The input "w" should be a function handle that
% accepts a scalar input and returns a scalar output. The default value
% of "w" is effectively w(x) = 1. The input "w" is specified using the
% WeightingFunction named argument.
%
% Example Usage:
%   [nodes, weights] = clenshawCurtis(N, a, b);
%   [nodes, weights] = clenshawCurtis(N, a, b, WeightingFunction=fun);
%   [nodes, weights] = clenshawCurtis(N, a, b, ComputeErrorWeights=false);
%   q = sum(fun(nodes) .* weights, 1);
%
% Inputs:
%   orderN - Scalar number of nodes to calculate.
%   a - Scalar integration lower bound. Must be real and finite.
%   b - Scalar integration upper bound. Must be real and finite.
% Outputs:
%   nodes - Column vector of coordinates at which to evaluate function.
%   weights - Column vector of weights to perform weighted sum.
%   errorWeights - Column vector of weights to estimate integration error.
% Named Options:
%   WeightingFunction - Optional weighting function w(x). See above.
%
% Author: Matt Dvorsky

arguments
   orderN(1, 1) {mustBeInteger, mustBePositive} = 10;
    a(1, 1) {mustBeReal, mustBeFinite} = -1;
    b(1, 1) {mustBeReal, mustBeFinite} = 1;
    options.WeightingFunction;
end

%% Calculate Moments
% Round orderN up to nearest even integer
N = 2*ceil(0.5*orderN);

n(:, 1) = 0:N;
if isfield(options, "WeightingFunction")
    moments = zeros(size(n));
    for ii = 1:length(moments)
        moments(ii) = integral(...
            @(x) options.WeightingFunction(0.5*(b - a) .* cos(x) + 0.5*(a + b)) ...
            .* sin(x) .* cos((ii - 1)*x), 0, pi);
    end
else
    moments = 2 ./ (1 - n.^2);
    moments(2:2:end) = 0;
end

%% Calculate Weights using DCT Type 1
% Use even ifft to calculate DCT-1.
weights = ifft([moments; moments(end - 1:-1:2)]);
if all(isreal(moments))
    weights = real(weights);
end
weights = weights(1:N + 1);
weights(2:end - 1) = 2*weights(2:end - 1);

%% Calculate Nodes
nodes(:, 1) = cos((0:N) * (pi/N));

%% Change Interval
weights = 0.5*(b - a) .* weights;
nodes = 0.5*(b - a) .* nodes + 0.5*(a + b);

end

