function [A1, A2, b1, b2] = constructMatrixEquation(O, nLayerInt)
%CONSTRUCTMATRIXEQUATION Computes the matrices A1, A2, b1, b2
% Computes the matrices A1, A2, b1, b2, used to solve the rectangular 
% waveguide equation for multilayer structures.
%
% Inputs:
%   nLayerInt - The values of the integrals I_ii(m, n, p, q). The size of
%       nLayerInt should be n-by-n-by-1-by-..., where n is the number of 
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
kc0m(1, :) = O.modeCutoffs;
kc0n(:, 1) = O.modeCutoffs;

%% Assemble Matrix Equation
% A1: Dependent on I(m, n) values (i.e., nLayerInt).
A1 = kc0m.^2 .* besselj(0, O.r * kc0m) .* nLayerInt(:, :, 1, :);

% A2: Independent of I(m, n) values. Defines exitation mode.
A2 = 0.5 * diag(besselj(0, O.r * kc0n));

% b1: Dependent on I(m, n) values. Always equal to first column of A1.
b1 = kc0m.^2 .* besselj(0, O.r * kc0m) .* nLayerInt(:, 1, 1, :);

% b2: Independent of II(m, n) values. Defines exitation mode
b2 = zeros(size(kc0n));
b2(1) = 0.5 * besselj(0, O.r * kc0n(1));

end


