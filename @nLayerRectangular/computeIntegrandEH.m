function [integrandE, integrandH] = computeIntegrandEH(O, tauP)
%computeIntegrandH Compute the partial integrand for I(n, q).
% This function can be used to compute I^(e)_ii(m, n, p, q) and 
% I^(h)_ii(m, n, p, q) by integrating the product of this function and
% the output of the "multilayerSpectrumEH" function over the interval 
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
% Note that the "multilayerSpectrumEH" function uses tau as its
% integration varable, while this function uses tauP. The following
% example shows the calculation of the integrals I^(e)_ii(...) and
% I^(h)_ii(...).
%
% Example:
%   function y = f(tauP)
%       tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
%       [specE, specH] = multilayerSpectrumEH(tau, k0, er, ur, thk);
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

%% Compute Integrals Over Psi at All Values of Tau
[psi(1, 1, 1, 1, :), weightsPsi(1, 1, 1, 1, :)] = ...
    O.fejer2(O.integralPointsPsi, 0, 0.5*pi);
xi = tau .* cos(psi);
eta = tau .* sin(psi);

fsc = 4*sum(weightsPsi .* tau.^3 .* sin(psi).^2 .* cos(psi).^2 .* f(xi, eta), 5);
fss = 4*sum(weightsPsi .* tau.^3 .* sin(psi).^4 .* f(xi, eta), 5);
fcc = 4*sum(weightsPsi .* tau.^3 .* cos(psi).^4 .* f(xi, eta), 5);

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
integrandE(tauP == 0, :, :, :) = 0;
integrandH(tauP == 0, :, :, :) = 0;

integrandE(tauP == 1, :, :, :) = 0;
integrandH(tauP == 1, :, :, :) = 0;

end
