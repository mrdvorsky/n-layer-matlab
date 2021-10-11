function [nodes, weights] = gaussLegendre(orderN, a, b)
%GAUSSLEGENDRE Generate Gauss-Legendre weights and nodes for closed interval integration.
%   This function generates the weights and nodes required to compute a
%   definite integral over a closed interval. The weights and nodes are
%   defined using the Gauss-Legendre Quadrature rules.
%   
%   The function outputs "nodes" and "weights" can be used to approximate
%   the definite integral of a function f(x)dx over the interval [a,b] by
%   computing I = sum(weights .* f(nodes)). This should give approximately
%   the same result as I = integral(f, a, b), with a higher value of
%   orderN resulting in a better approximation. The parameter orderN is
%   the number of points at which to evaluate f(x). If f(x) is a polynomial
%   with degree less than to 2*orderN, the result will be exact.
%   The inputs "a" and "b" must be real, finite values.

%% Check Input
if nargin < 3
    b = 1;
end
if nargin < 2
    a = -1;
end
if nargin < 1
    orderN = 10;
end

%% Calculate Weights and Nodes
N=orderN-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
nodes=(-1*(1-y)+1*(1+y))/2;      

% Compute the weights
weights=(1-(-1))./((1-y.^2).*Lp.^2)*(N2/N1)^2;

%% Change Interval
weights = 0.5*(b - a) .* weights;
nodes = 0.5*(b - a) .* nodes + 0.5*(a + b);

end

