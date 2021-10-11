function [A1, A2, b1, b2] = constructMatrixEquation(O, nLayerInt)
%CONSTRUCTMATRIXEQUATION Summary of this function goes here
%   Detailed explanation goes here
%
%   Designation for each dimension:
%       1: Mode matrix columns
%       2: Mode matrix rows
%       3: Subscript i of integration parameters I_i(n,m,p,q)

%% Mode Coefficients
am(1, :) = O.modes(:, 1) * pi ./ O.a;
bn(1, :) = O.modes(:, 2) * pi ./ O.b;
ap(:, 1) = O.modes(:, 1) * pi ./ O.a;
bq(:, 1) = O.modes(:, 2) * pi ./ O.b;

indTE = find(O.modes(:, 1) > 0);
indTM = find(O.modes(:, 2) > 0);

%% Assemble Matrix Equation
% A1: Dependent on I(m,n,p,q) values (i.e., nLayerInt)
A1_EE = (bn(1, indTE) .* nLayerInt(indTE, indTE, 1, :) ...
    - am(1, indTE) .* nLayerInt(indTE, indTE, 2, :));
A1_EM = (am(1, indTM) .* nLayerInt(indTE, indTM, 1, :) ...
    + bn(1, indTM) .* nLayerInt(indTE, indTM, 2, :));
A1_ME = (bn(1, indTE) .* nLayerInt(indTM, indTE, 3, :) ...
    - am(1, indTE) .* nLayerInt(indTM, indTE, 4, :));
A1_MM = (am(1, indTM) .* nLayerInt(indTM, indTM, 3, :) ...
    + bn(1, indTM) .* nLayerInt(indTM, indTM, 4, :));

A1 = [A1_EE, A1_EM; ...
    A1_ME, A1_MM];

% A2: Independent of I(m,n,p,q) values
A2_EE = -diag(ap(indTE, 1) .* (1 + (bq(indTE, 1) == 0)));
A2_MM = diag(ap(indTM, 1));

A2_EM = zeros(length(indTE), length(indTM));
A2_ME = zeros(length(indTM), length(indTE));
A2_EM(indTM, :) = diag(bq(indTM, 1));
A2_ME(:, indTM) = diag(bq(indTM, 1) .* (1 + (ap(indTM, 1) == 0)));

A2 = [A2_EE, A2_EM; ...
    A2_ME, A2_MM] .* (0.25 * O.a * O.b);

% b1: Dependent on I(m,n,p,q) values
b1_E = ap(1) .* nLayerInt(indTE, 1, 2, :);
b1_M = ap(1) .* nLayerInt(indTM, 1, 4, :);

b1 = [b1_E; b1_M];

% b2: Independent of I(m,n,p,q) values
b2 = zeros(length(indTE) + length(indTM), 1);
b2(1) = b2(1) - 2*ap(1) .* (0.25 * O.a * O.b);

end


