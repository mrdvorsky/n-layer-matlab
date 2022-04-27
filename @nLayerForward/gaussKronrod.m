function [nodes, weights, errorWeights] = gaussKronrod(numSegs, a, b)
%GAUSSLEGENDRE Generate Gauss-Kronrod weights and nodes for closed interval integration.
% This function generates the weights and nodes required to compute a
% definite integral over a closed interval. The weights and nodes are
% defined using the Gauss-Kronrod Quadrature rules.
%
% The function outputs "nodes" and "weights" can be used to approximate
% the definite integral of a function f(x)dx over the interval [a,b] by
% computing q = sum(weights .* f(nodes)). This should give approximately
% the same result as q = integral(f, a, b), with a higher value of
% numSegs resulting in a better approximation. The output "errorWeights"
% can be similarly used to estimate the error in the evaluation of q (i.e,
% q = sum(errorWeights .* f(nodes)); ).
%
% Use of these quadrature rules will never result in evalation of the
% function at the interval endpoints a and b.
%
% Example Usage:
%   [nodes, weights, errorWeights] = gaussKronrod(N, a, b);
%   q = sum(fun(nodes) .* weights, 1);
%   qErr = sum(fun(nodes) .* errorWeights, 1);
%
% Inputs:
%   numSegs - Number of uniformly-sized segments to subdivide [a, b].
%   a - Scalar integration lower bound. Must be real and finite.
%   b - Scalar integration upper bound. Must be real and finite.
% Outputs:
%   nodes - Column vector of coordinates at which to evaluate function.
%   weights - Column vector of weights to perform weighted sum.
%   errorWeights - Column vector of weights to perform weighted sum for
%       estimation of error.
%
% Author: Matt Dvorsky

arguments
    numSegs(1, 1) {mustBeInteger, mustBePositive} = 10;
    a(1, 1) {mustBeReal, mustBeFinite} = -1;
    b(1, 1) {mustBeReal, mustBeFinite} = 1;
end

%% Generate Gauss-Kronrod Weights and Nodes
% Gauss-Kronrod (7,15) pair. Use symmetry in defining nodes and weights.
pnodes = [ ...
    0.2077849550078985; 0.4058451513773972; 0.5860872354676911; ...
    0.7415311855993944; 0.8648644233597691; 0.9491079123427585; ...
    0.9914553711208126];
pwt = [ ...
    0.2044329400752989; 0.1903505780647854; 0.1690047266392679; ...
    0.1406532597155259; 0.1047900103222502; 0.06309209262997855; ...
    0.02293532201052922];
pwt7 = [0; 0.3818300505051189; 0; ...
    0.2797053914892767; 0; 0.1294849661688697; 0];

nodesGK = [-pnodes(end:-1:1); 0; pnodes];
weightsGK = [pwt(end:-1:1); 0.2094821410847278; pwt];
errorWeightsGK = weightsGK - [pwt7(end:-1:1); 0.4179591836734694; pwt7];

%% Split Interval and Assign Weights and Nodes
intervalPoints = linspace(a, b, numSegs + 1);
intervals = [intervalPoints(1:end - 1); intervalPoints(2:end)];

midpoints = 0.5 * sum(intervals, 1);
halfLengths = 0.5 * diff(intervals, 1);

nodes = reshape(nodesGK .* halfLengths + midpoints, [], 1);
weights = reshape(weightsGK .* halfLengths, [], 1);
errorWeights = reshape(errorWeightsGK .* halfLengths, [], 1);

end

