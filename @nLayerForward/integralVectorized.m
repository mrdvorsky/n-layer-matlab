function [q] = integralVectorized(fun, a, b, options)
%INTEGRALVECTORIZED Numerically evaluate integral using adaptive Gauss-Kronrod quadrature.
%Based on quadgk but modified to work quickly for array-valued functions.
% Example Usage:
%   q = integralVectorized(fun, a, b);
%   q = integralVectorized(fun, a, b, RelTol=1e-4);
%   q = integralVectorized(fun, a, b, InitialIntervalCount=20);
%   q = integralVectorized(fun, a, b, MaxFunctionEvaluations=20000);
%
% q = integralVectorized(fun, a, b) attempts to approximate the integral 
% of an array-valued function FUN from A to B using high order global 
% adaptive quadrature and default error tolerances.
%
% The function Y = FUN(X) should accept a column vector argument X and 
% return an array result Y, the columns of which are the integrand 
% evaluated at each element of Y. The size of the first dimension of Y 
% should match X, the remaining dimensions can be any size.
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
% The size of the output Q will be the same as FUN((A + B)/2).
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
%   ErrInd (1): Indices ii of q(:, ii) to use for error bounds. See above.
%
% Author: Matt Dvorsky

arguments
    fun;
    a(1, 1) {mustBeNumeric, mustBeFinite} = -1;
    b(1, 1) {mustBeNumeric, mustBeFinite} = 1;
    
    options.RelTol(1, 1) {mustBeNumeric, mustBeFinite} = 1e-6;
    options.AbsTol(1, 1) {mustBeNumeric, mustBeFinite} = 1e-8;
    options.Verbosity double {mustBeNumeric, mustBeFinite} = 0;
    options.InitialIntervalCount {mustBeNumeric, mustBeFinite} = 9;
    options.MaxIntervalCount {mustBeNumeric} = 2000;
    options.MaxFunctionEvaluations {mustBeNumeric} = 30000;
    options.ErrInd(:, 1) = 1;
end

%% Generate Gauss-Kronrod Weights and Nodes
% Gauss-Kronrod (7,15) pair. Use symmetry in defining nodes and weights.
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
qSize = 0;

numEvaluations = 0;

%% Main Integration Loop
while true
    % Calculate new function evaluation coordinates
    midpoints = 0.5 * sum(intervals, 1);
    halfLengths = 0.5 * diff(intervals, 1);
    x = reshape(nodesGK .* halfLengths + midpoints, [], 1);
            
    % Evaluate the function at new coordinates
    fx = fun(x);
    qSize = size(fx);
    fx = reshape(fx, length(weightsGK), length(midpoints), []);
    numEvaluations = numEvaluations + length(x);
    
    % Compute integral and error estimate for each subinterval
    qIntervals = sum(weightsGK .* halfLengths .* fx, 1);
    errorIntervals = sum(errorWeightsGK .* halfLengths .* fx(:, :, options.ErrInd), 1);
    
    % Calculate current values of q
    q = qPartial + sum(qIntervals, 2);
    
    % Calculate absolute tolerance from relative tolerance
    AbsTolTmp = max(options.RelTol .* abs(q(1, 1, options.ErrInd)), options.AbsTol);
    
    % Find indices of the subinterals that are sufficiently accurate
    accurateIntervals = all(abs(errorIntervals) <= ...
        (2*AbsTolTmp ./ pathlength) .* abs(halfLengths), 3);
    
    % Update cumulative error of all sufficiently accurate intervals
    errorPartial = errorPartial + sum(errorIntervals(1, accurateIntervals, :), 2);
    
    % Calculate upper bound of current error estimate
    errorBound = abs(errorPartial) ...
        + sum(abs(errorIntervals(1, ~accurateIntervals, :)), 2);
    
    % Error if we have non-finite values
    if ~(all(isfinite(q)) && all(isfinite(errorBound)))
        error("Calculation returned a non-finite value.");
    end
    
    % Check for convergence or if all intervals are accurate
    if all(errorBound <= AbsTolTmp) || all(accurateIntervals)
        if options.Verbosity > 0
            fprintf("Converged after (%d) function evaluations.\n", ...
                numEvaluations);
        end
        break;
    end
    
    % Update cumulative integral of all sufficiently accurate intervals
    qPartial = qPartial + sum(qIntervals(1, accurateIntervals, :), 2);
    
    % Split the remaining subintervals in half
    intervals = reshape( ...
        [intervals(1, ~accurateIntervals); ...
        midpoints(~accurateIntervals); ...
        midpoints(~accurateIntervals); ...
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

%% Reshape Output
q = reshape(q, [1, qSize(2:end)]);

end

