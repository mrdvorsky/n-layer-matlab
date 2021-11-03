function [integrandH] = computeIntegrandH(O, tauP)
%computeIntegrandEH Compute the partial integrand for I(m, n, p, q).
% This function can be used to compute I^(e)_ii(m, n, p, q) and 
% I^(h)_ii(m, n, p, q) by integrating the product of this function and
% the output of the "multilayerSpectrumRect" function over the interval 
% from 0 to 1. This integral can then be used with the
% "constructMatrixEquation" function to construct the matrix equation 
% required to find to solve for the reflection coefficient. See
% documentation of "constructMatrixEquation" for more details.
%
% This function computes the integrand for I^(e)_ii(...) and
% I^(h)_ii(...) after a change of variables from tau to tauP.
% The variable tau is related to tauP by the following equation.
%       tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
%
% Note that the "multilayerSpectrumRect" function uses tau as its
% integration varable, while this function uses tauP. The following
% example shows the calculation of the integrals I^(e)_ii(...) and
% I^(h)_ii(...).
%
% Example:
%   function y = f(tauP)
%       tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
%       [specE, specH] = multilayerSpectrumRect(tau, k0, er, ur, thk);
%       [integrandE, integrandH] = O.computeIntegrandEH(tauP);
%       y = specE.*integrandE + specH.*integrandH;
%   end
%   nLayerInt = integral(f, 0, 1, "ArrayValued", "true");
%
% The output "nLayerInt" can then be used as an input to the
% "constructMatrixEquation" function. See documentation for more details.
%
% Designation for each dimension:
%   1: Integration variable tauP
%   2: Mode matrix columns
%   3: Mode matrix rows
%   4: Subscript i of integration parameters I_ii(m, n, p, q)
%   5: Integration variable psi (temporarily)
%
% Author: Matt Dvorsky

arguments
    O;
    tauP(:, 1);
end

%% Compute Interpolating Coordinates
% The integral needs to be evaulated from tau = [0, inf). However, a change
% of variables tau = L(1 - tauP)/tauP is used here so that the interpolant
% can be uniform in (0, 1].
L = O.integralScaleFactor;
tau = L * (1 - tauP) ./ tauP;

% Weighting function to account for change of variables.
weights = L ./ (tauP.^2);

%% Compute Waveguide Mode Cutoffs
kc0m(1, 1, :, 1) = O.modeCutoffs;
kc0n(1, :, 1, 1) = O.modeCutoffs;

%% Compute Integrand Values
integrandH = (weights .* tau) .* (besselj(1, O.r * tau).^2 ./ ...
    ((tau.^2 - (kc0n).^2) .* (tau.^2 - (kc0m).^2)));

%% Fix Nans Caused by Singularities
integrandH(tauP == 0, :, :, :) = 0;
integrandH(tauP == 1, :, :, :) = 0;

end
