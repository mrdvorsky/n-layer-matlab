function [nodes, weights] = gaussLegendre(orderN, a, b)
%GAUSSLEGENDRE Generate Gauss-Legendre weights and nodes for closed interval integration.
% This function generates the weights and nodes required to compute a
% definite integral over a closed interval. The weights and nodes are
% defined using the Gauss-Legendre Quadrature rules.
%
% The function outputs "nodes" and "weights" can be used to approximate
% the definite integral of a function f(x)dx over the interval [a,b] by
% computing q = sum(weights .* f(nodes)). This should give approximately
% the same result as q = integral(f, a, b), with a higher value of
% orderN resulting in a better approximation. The parameter orderN is
% the number of points at which to evaluate f(x).
%
% If f(x) is a polynomial with degree less than to 2*orderN, the result
% will be exact.
%
% Use of these quadrature rules will never result in evalation of the
% function at the interval endpoints a and b.
%
% The weights and nodes are computed in two steps. First, the Golub-Welsh
% algorithm is used to compute the nodes. This step involves finding the
% eigenvalues of a symmetric tridiagonal matrix. Typically, the
% eigenvectors of this matrix are used to compute the wieghts; however,
% this requires O(n^2) memory, which is prohibitive. Instead, only the
% eigenvalues are computed using O(n) memory. The next step computes the
% weights by using the three-term recursion formula for the Legendre
% polynomials. Each step requires O(n) memory and O(n^2) time.
%
% Example Usage:
%   [nodes, weights] = gaussLegendre(N, a, b);
%   q = sum(fun(nodes) .* weights, 1);
%
% Inputs:
%   orderN - Scalar number of nodes to calculate.
%   a - Scalar integration lower bound. Must be finite.
%   b - Scalar integration upper bound. Must be finite.
% Outputs:
%   nodes - Column vector of coordinates at which to evaluate function.
%   weights - Column vector of weights to perform weighted sum.
%
% Author: Matt Dvorsky

arguments
    orderN(1, 1) {mustBeInteger} = 10;
    a(1, 1) {mustBeNumeric} = -1;
    b(1, 1) {mustBeNumeric} = -1;
end

%% Calculate Nodes
% Golub-Welsh Algorithm
alpha(:, 1) = 1:(orderN - 1);
beta = alpha ./ sqrt(4*alpha.^2 - 1);
J = spdiags([0; beta], 1, orderN, orderN) + spdiags([beta; 0], -1, orderN, orderN);

% The nodes are the eigenvalues of the the matrix J
nodes = eig(J);

%% Calculate Weights Using Legendre Polynomial
% The weights can be computed by evaluating the derivative of P_n(x) at
% each node coordinate, where P_n(x) is the nth order Legendre polynomial.
% P_n'(x) can be evaluated using the recurrence relationship.

% Initial Values of P_0(x), P_0'(x), P_1(x), P_1'(x)
pn_prev = 1 + 0*nodes;      % P_0(x) = 1
pn = nodes;                 % P_1(x) = x

pn_prime_prev = 0*nodes;    % P_0'(x) = 0
pn_prime = 1 + 0*nodes;     % P_1'(x) = 1
for n = 2:orderN
    pn_next = ((2*n - 1).*nodes.*pn - (n - 1).*pn_prev) ./ n;
    pn_deriv_next = ((2*n - 1).*(pn + nodes.*pn_prime) ...
        - (n - 1).*pn_prime_prev) ./ n;
    
    % Update For Next Loop
    [pn_prev, pn] = deal(pn, pn_next);
    [pn_prime_prev, pn_prime] = deal(pn_prime, pn_deriv_next);
end

% Compute weights from P_n'(x)
weights = 2 ./ (1 - nodes.^2) ./ pn_prime.^2;

%% Change Interval
weights = 0.5*(b - a) .* weights;
nodes = 0.5*(b - a) .* nodes + 0.5*(a + b);

end

