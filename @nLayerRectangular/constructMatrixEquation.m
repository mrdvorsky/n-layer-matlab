function [A1, A2, b1, b2] = constructMatrixEquation(O, nLayerInt)
%CONSTRUCTMATRIXEQUATION Computes the matrices A1, A2, b1, b2
% Computes the matrices A1, A2, b1, b2, used to solve the rectangular 
% waveguide equation for multilayer structures.
%
% Inputs:
%   nLayerInt - The values of the integrals I_ii(m, n, p, q). The size of
%       nLayerInt should be n-by-n-by-4-by-..., where n is the number of 
%       TE modes to consider. The value of nLayerInt(kk, ll, ii) should be 
%       equal to I_ii(m, n, p, q), where [m, n] = modesTE(kk, :), and 
%       [p, q] = modesTE(ll, :). The 4th dimensions onward can be any size
%       and will be preserved in the outputs.
% Outputs:
%   A1, A2, b1, b2 - Matrices used to compute S11. See usage below. The
%       4th dimensions and onward of the outputs will match the dimensions
%       of the input nLayerInt.
%
% The outputs of this function can be used along with the outputs of the
% constructFrequencyMultipliers(...) function to calculate S11 for a
% rectangular waveguide. See example usage below.
%
% Example Usage 1:
%   [k_A1, k_A2, k_b1, k_b2] = O.constructFrequencyMultipliers(f);
%   [A1, A2, b1, b2] = O.constructMatrixEquation(nLayerInt);
%   sqrt(ur(1) ./ er(1));
%   x = (A1.*k_A1 + etaR1.*A2.*k_A2) ...
%        \ (b1.*k_b1 + etaR1.*b2.*k_b2);
%   S11 = x(1);
%
% Alternatively, this function can be used with no arguments to compute
% only A2 and b2, and A1 and b1 can be computed using the "computeA1b1"
% function.
%
% Example Usage 2:
%   ...
%   [~, A2, ~, b2] = O.constructMatrixEquation();
%   [A1, b1] = computeA1b1(f, er, ur, thk, AbsTol);
%   ...
%
% The input parameter nLayerInt can be calculated by making use of the
% "computeIntegrandEH" function. See documentation for
% "computeIntegrandEH" for more details.
%
% Designation for each dimension (4th dimension and higher are duplicated):
%   1: Mode matrix columns
%   2: Mode matrix rows
%   3: Subscript ii of integration parameters I_ii(m, n, p, q)
%
% Author: Matt Dvorsky

arguments
    O;
    nLayerInt = [];
end

%% Mode Coefficients
am(1, :) = O.modesTE(:, 1) * pi ./ O.a;
bn(1, :) = O.modesTE(:, 2) * pi ./ O.b;
ap(:, 1) = O.modesTE(:, 1) * pi ./ O.a;
bq(:, 1) = O.modesTE(:, 2) * pi ./ O.b;

indTE = find(O.modesTE(:, 1) > 0);
indTM = find(O.modesTE(:, 2) > 0);

%% Assemble Matrix Equation
% A1: Dependent on I(m, n, p, q) values (i.e., nLayerInt)
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

% A2: Independent of I(m, n, p, q) values. Defines exitation mode
A2_EE = -diag(ap(indTE, 1) .* (1 + (bq(indTE, 1) == 0)));
A2_MM = diag(ap(indTM, 1));

A2_EM = zeros(length(indTE), length(indTM));
A2_ME = zeros(length(indTM), length(indTE));
A2_EM(indTM, :) = diag(bq(indTM, 1));
A2_ME(:, indTM) = diag(bq(indTM, 1) .* (1 + (ap(indTM, 1) == 0)));

A2 = [A2_EE, A2_EM; ...
    A2_ME, A2_MM] .* (0.25 * O.a * O.b);

% b1: Dependent on I(m, n, p, q) values
b1_E = ap(1) .* nLayerInt(indTE, 1, 2, :);
b1_M = ap(1) .* nLayerInt(indTM, 1, 4, :);

b1 = [b1_E; b1_M];

% b2: Independent of I(m, n, p, q) values. Defines exitation mode
b2 = zeros(length(indTE) + length(indTM), 1);
b2(1) = b2(1) - 2*ap(1) .* (0.25 * O.a * O.b);

end


