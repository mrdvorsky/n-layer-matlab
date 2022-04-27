function [A1, A2] = constructMatrixEquation(O, nLayerInt)
%CONSTRUCTMATRIXEQUATION Computes the matrices A1 and A2
% Computes the matrices A1 and A2 used to solve the circular TE0n mode 
% waveguide equation for multilayer structures.
%
% Inputs:
%   nLayerInt - The values of the integrals I(m, n). The size of nLayerInt
%       should be n-by-n-by-..., where n is the number of TE0n modes to
%       consider. The value of nLayerInt(kk, ll, :) should be equal to I(m,
%       n), where m and n are elements of modesTE. The 3rd dimensions
%       onward can be any size and will be preserved in the outputs.
% Outputs:
%   A1 and A2 - Matrices used to compute S11. See usage below. The 3rd
%       dimensions and onward of the outputs will match the dimensions
%       of the input nLayerInt.
%
% The outputs of this function can be used along with the outputs of the
% constructFrequencyMultipliers(...) function to calculate S11 for a
% rectangular waveguide. See example usage below.
%
% Example Usage 1:
%   [k_A2, k_b2] = O.constructFrequencyMultipliers(f);
%   [A1, A2] = O.constructMatrixEquation(nLayerInt);
%   x = (A1 + ur.*A2.*k_A2) \ (-A1(:, 1) + ur.*A2(:, 1).*k_b2);
%   S11 = x(1);
%
% Alternatively, this function can be used with no arguments to compute
% only A2, and A1 can be computed using the "computeA1b1" function.
%
% Example Usage 2:
%   ...
%   [~, A2] = O.constructMatrixEquation();
%   A1 = computeA1(f, er, ur, thk);
%   ...
%
% The input parameter nLayerInt can be calculated by making use of the
% "computeIntegrandH" function. See documentation for "computeIntegrandH"
% for more details.
%
% Designation for each dimension (3rd dimension and higher are duplicated):
%   1: Mode matrix columns
%   2: Mode matrix rows
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
A1 = kc0m.^2 .* besselj(0, O.r * kc0m) .* nLayerInt;

% A2: Independent of I(m, n) values. Defines exitation mode.
A2 = 0.5 * diag(besselj(0, O.r * kc0n));

end


