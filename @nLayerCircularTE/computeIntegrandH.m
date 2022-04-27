function [integrandH] = computeIntegrandH(O, tauP)
%computeIntegrandH Compute the partial integrand for I(m, n).
% This function can be used to compute I(m, n) by integrating the product
% of this function and the output of the "multilayerSpectrumH" function
% over the interval from 0 to 1. This integral can then be used with the
% "constructMatrixEquation" function to construct the matrix equation 
% required to find to solve for the reflection coefficient. See
% documentation of "constructMatrixEquation" for more details.
%
% This function computes the integrand for I(...) after a change of
% variables from tau to tauP. The variable tau is related to tauP by the
% following equation.
%       tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
%
% Note that the "multilayerSpectrumH" function uses tau as its
% integration varable, while this function uses tauP. The following
% example shows the calculation of the integral I(...).
%
% Example:
%   function y = f(tauP)
%       tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
%       specH = multilayerSpectrumH(tau, k0, er, ur, thk);
%       integrandH = O.computeIntegrandH(tauP);
%       y = specH.*integrandH;
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
%
% Author: Matt Dvorsky

arguments
    O;
    tauP(:, 1);
end

%% Compute Interpolating Coordinates
% The integral needs to be evaluated from tau = [0, inf). However, a change
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
integrandH(tauP == 0, :, :) = 0;
integrandH(tauP == 1, :, :) = 0;

end
