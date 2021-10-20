function [q, nodes_out, weights_out] = integralWeightsAndNodes(fun, a, b, options)
%integralWeightsAndNodes Numerically evaluate integral using adaptive Gauss-Kronrod quadrature.
%Return final weights and nodes used as well as the computed integral.
% Example Usage:
%   [q, nodes, weights] = integralWeightsAndNodes(fun, a, b);
%   [q, nodes, weights] = integralWeightsAndNodes(fun, a, b, RelTol=1e-4);
%
% q = integralVectorized(fun, a, b) attempts to approximate the integral 
% of an array-valued function FUN from A to B using high order global 
% adaptive quadrature and default error tolerances.
%
% The function Y = FUN(X) should accept a column vector argument X and 
% return an column vector result Y with the same size as X evaluated at 
% each element of Y.
%
% FUN must be a function handle. The parameters A and B
% must both be finite and real and in ascending order.
%
% This function attempts to satisfy the following relationship:
%       ERRBND(errInd) <= max(RELTOL*|Q(errInd)|, ABSTOL)
% This constraint is the same as the one used in the built-in "integral"
% funtion. The default value of RELTOL is 1e-6, and the default of ABSTOL
% is 1e-8.
%
% The output parameters nodes_out and weights_out are the final nodes and
% weights used. In other words, q = sum(fun(nodes) .* weights);
%
% Optional Arguments:
%   RelTol (1e-6): Relative function tolerance constraint. See above.
%   AbsTol (1e-8): Absolute function tolerance constraint. See above.
%   Verbosity (0): Set to 1 or higher for console output upon convergence.
%   InitialIntervalCount (2000): Maximum number of subintervals to attempt
%       before giving up.
%   MaxIntervalCount (10000): Maximum number of subintervals to consider
%       simultaneously before giving up.
%   MaxFunctionEvaluations (30000): Maximum number of function evaluations
%       to attempt before giving up.
%
% Author: Matt Dvorsky

arguments
    fun;
    a(1, 1) {mustBeNumeric, mustBeFinite} = -1;
    b(1, 1) {mustBeNumeric, mustBeFinite} = 1;
    
    options.RelTol(1, 1) {mustBeNumeric, mustBeFinite} = 1e-6;
    options.AbsTol(1, 1) {mustBeNumeric, mustBeFinite} = 1e-8;
    options.Verbosity {mustBeNumeric, mustBeFinite} = 0;
    options.InitialIntervalCount {mustBeNumeric, mustBeFinite} = 9;
    options.MaxIntervalCount {mustBeNumeric} = 2000;
    options.MaxFunctionEvaluations {mustBeNumeric} = 30000;
end

%% Generate Gauss-Kronrod Weights and Nodes
% Gauss-Kronrod (7,15) pair.
nodesGK = [ ...
    -0.9914553711208126; -0.9491079123427585; -0.8648644233597691; ...
    -0.7415311855993944; -0.5860872354676911; -0.4058451513773972; ...
    -0.2077849550078985;  0;                   0.2077849550078985; ...
     0.4058451513773972;  0.5860872354676911;  0.7415311855993944; ...
     0.8648644233597691;  0.9491079123427585;  0.9914553711208126];
weightsGK = [ ...
    0.0229353220105292; 0.0630920926299786; 0.1047900103222502; ...
    0.1406532597155259; 0.1690047266392679; 0.1903505780647854; ...
    0.2044329400752989; 0.2094821410847278; 0.2044329400752989; ...
    0.1903505780647854; 0.1690047266392679; 0.1406532597155259; ...
    0.1047900103222502; 0.0630920926299786; 0.0229353220105292];
errorWeightsGK = [ ...
     0.0229353220105292; -0.0663928735388910;  0.1047900103222502; ...
    -0.139052131773751;   0.1690047266392679; -0.191479472440334; ...
     0.2044329400752989; -0.2084770425887420;  0.2044329400752989; ...
    -0.191479472440334;   0.1690047266392679; -0.139052131773751; ...
     0.1047900103222502; -0.0663928735388910;  0.0229353220105292];

%% Initialize Integral Loop
% INTERVALS contains subintervals of [a,b] where the integral is not
%  accurate enough. First and second rows are left and right endpoints.
intervalPoints = linspace(a, b, options.InitialIntervalCount + 1);
intervals = [intervalPoints(1:end - 1); intervalPoints(2:end)];

pathlength = b - a;
qPartial = 0;
errorPartial = 0;
numEvaluations = 0;

nodes_out = [];
weights_out = [];

%% Main Integration Loop
while true
    % Calculate new function evaluation coordinates
    intervalMidpoints = 0.5 * sum(intervals, 1);
    intervalsHalfLength = 0.5 * diff(intervals(:, 1), 1);
    weights = weightsGK .* intervalsHalfLength;
    errorWeights = errorWeightsGK .* intervalsHalfLength;
    x = reshape(nodesGK .* intervalsHalfLength + intervalMidpoints, ...
        length(weights), length(intervalMidpoints));
    
    % Evaluate the function at new coordinates
    fx = reshape(fun(x(:)), length(weights), length(intervalMidpoints));
    numEvaluations = numEvaluations + numel(x);
    
    % Compute integral and error estimate for each subinterval
    qIntervals = sum(weights .* fx, 1);
    errorIntervals = sum(errorWeights .* fx, 1);
    
    % Calculate current values of q
    q = qPartial + sum(qIntervals, 2);
    
    % Calculate absolute tolerance from relative tolerance
    AbsTolTmp = max(options.RelTol .* abs(q), options.AbsTol);
    
    % Find indices of the subinterals that are sufficiently accurate
    accurateIntervals = all(abs(errorIntervals) <= ...
        (2*AbsTolTmp ./ pathlength) .* abs(intervalsHalfLength), 3);
    
    % Update cumulative error of all sufficiently accurate intervals
    errorPartial = errorPartial + sum(errorIntervals(accurateIntervals));
    
    % Calculate upper bound of current error estimate
    errorBound = abs(errorPartial) ...
        + sum(abs(errorIntervals(~accurateIntervals)));
    
    % Error if we have non-finite values
    if ~(isfinite(q) && isfinite(errorBound))
        error("Calculation returned a non-finite value.");
    end
    
    % Check for convergence
    if errorBound <= AbsTolTmp
        accurateIntervals(:) = true;
    end
    
    nodes_out = [nodes_out; reshape(x(:, accurateIntervals), [], 1)];
    weights_out = [weights_out; repmat(weights, sum(accurateIntervals), 1)];
    
    if all(accurateIntervals)
        if options.Verbosity > 0
            fprintf("Converged after (%d) function evaluations.\n", ...
                numEvaluations);
        end
        break;
    end
    
    % Update cumulative integral of all sufficiently accurate intervals
    qPartial = qPartial + sum(qIntervals(accurateIntervals));
    
    % Split the remaining subintervals in half
    intervals = reshape( ...
        [intervals(1, ~accurateIntervals); ...
        intervalMidpoints(~accurateIntervals); ...
        intervalMidpoints(~accurateIntervals); ...
        intervals(2, ~accurateIntervals)], ...
        2, []);
    
    % Error if splitting results in too many subintervals
    if size(intervals, 2) > options.MaxIntervalCount
        error("Maximum number of subintervals reached (%d > %d).", ...
            size(intervals, 2), options.MaxIntervalCount);
    end
    
    % Error if splitting results in too many function evaluations
    if numel(nodesGK)*size(intervals, 2) + numEvaluations > options.MaxFunctionEvaluations
        error("Maximum number of subintervals reached (%d > %d).", ...
            size(intervals, 2), options.MaxIntervalCount);
    end
end

end




