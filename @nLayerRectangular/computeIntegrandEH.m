function [integrandE, integrandH] = computeIntegrandEH(O, kRhoP)
%computeIntegrandEH Compute the partial integrand for I(m, n, p, q).
% This function can be used to compute I^(e)_ii(m, n, p, q) and 
% I^(h)_ii(m, n, p, q) by integrating the product of this function and
% the output of the "multilayerSpectrumEH" function over the interval 
% from 0 to 1. This integral can then be used with the
% "constructMatrixEquation" function to construct the matrix equation 
% required to find to solve for the reflection coefficient. See
% documentation of "constructMatrixEquation" for more details.
%
% This function computes the integrand for I^(e)_ii(...) and
% I^(h)_ii(...) after a change of variables from kRho to kRhoP.
% The variable kRho is related to kRhoP by the following equation.
%       kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
%
% Note that the "multilayerSpectrumEH" function uses kRho as its
% integration varable, while this function uses kRhoP. The following
% example shows the calculation of the integrals I^(e)_ii(...) and
% I^(h)_ii(...).
%
% Example:
%   function y = f(kRhoP)
%       kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
%       [specE, specH] = multilayerSpectrumEH(kRho, k0, er, ur, thk);
%       [integrandE, integrandH] = O.computeIntegrandEH(kRhoP);
%       y = specE.*integrandE + specH.*integrandH;
%   end
%   nLayerInt = integral(f, 0, 1, "ArrayValued", "true");
%
% The output "nLayerInt" can then be used as an input to the
% "constructMatrixEquation" function. See documentation for more details.
%
% Designation for each dimension:
%   1: Integration variable kRhoP
%   2: Mode matrix columns
%   3: Mode matrix rows
%   4: Subscript i of integration parameters I_ii(m, n, p, q)
%   5: Integration variable kPhi (temporarily)
%
% Author: Matt Dvorsky

arguments
    O;
    kRhoP(:, 1);
end

%% Compute Interpolating Coordinates
% The integral needs to be evaluated from kRho = [0, inf). However, a change
% of variables kRho = L(1 - kRhoP)/kRhoP is used here so that the interpolant
% can be uniform in (0, 1].
L = O.integralScaleFactor;
kRho = L * (1 - kRhoP) ./ kRhoP;

% Weighting function to account for change of variables.
weights = L ./ (kRhoP.^2);

%% Compute Waveguide Mode Cutoffs
am(1, 1, :, 1) = O.modesTE(:, 1) * pi ./ O.a;
bn(1, 1, :, 1) = O.modesTE(:, 2) * pi ./ O.b;
ap(1, :, 1, 1) = O.modesTE(:, 1) * pi ./ O.a;
bq(1, :, 1, 1) = O.modesTE(:, 2) * pi ./ O.b;

%% Helper Functions
f = @(xi, eta) (O.a.^2 .* O.b.^2 / 16) ...
    .* sinc((0.5/pi)*O.a * (xi - am)) ...
    .* sinc((0.5/pi)*O.a * (xi - ap)) ...
    .* sinc((0.5/pi)*O.b * (eta - bn)) ...
    .* sinc((0.5/pi)*O.b * (eta - bq)) ...
    ./ ((am + xi) .* (ap + xi) ...
    .* (bn + eta) .* (bq + eta));

%% Compute Integrals Over kPhi at All Values of kRho
[kPhi(1, 1, 1, 1, :), weights_kPhi(1, 1, 1, 1, :)] = ...
    O.fejer2(O.integralPoints_kPhi, 0, 0.5*pi);
xi = kRho .* cos(kPhi);
eta = kRho .* sin(kPhi);

fsc = 4*sum(weights_kPhi .* kRho.^3 .* sin(kPhi).^2 .* cos(kPhi).^2 .* f(xi, eta), 5);
fss = 4*sum(weights_kPhi .* kRho.^3 .* sin(kPhi).^4 .* f(xi, eta), 5);
fcc = 4*sum(weights_kPhi .* kRho.^3 .* cos(kPhi).^4 .* f(xi, eta), 5);

%% Compute Integrand Values
integrandE = zeros(size(fsc)) + zeros(1, 1, 1, 4);
integrandH = zeros(size(integrandE));

integrandE(:, :, :, 1) = weights .* (bn .* ap .* (4/pi.^2) .* fsc);
integrandE(:, :, :, 2) = weights .* (am .* ap .* (4/pi.^2) .* fss);
integrandE(:, :, :, 3) = weights .* (bn .* bq .* (4/pi.^2) .* fcc);
integrandE(:, :, :, 4) = weights .* (am .* bq .* (4/pi.^2) .* fsc);

integrandH(:, :, :, 1) = -weights .* (bn .* ap .* (4/pi.^2) .* fsc);
integrandH(:, :, :, 2) =  weights .* (am .* ap .* (4/pi.^2) .* fsc);
integrandH(:, :, :, 3) =  weights .* (bn .* bq .* (4/pi.^2) .* fsc);
integrandH(:, :, :, 4) = -weights .* (am .* bq .* (4/pi.^2) .* fsc);

%% Fix Nans Caused by Singularities
integrandE(kRhoP == 0, :, :, :) = 0;
integrandH(kRhoP == 0, :, :, :) = 0;

integrandE(kRhoP == 1, :, :, :) = 0;
integrandH(kRhoP == 1, :, :, :) = 0;

end
