function [integralE, integralH] = computeIntegralEH(O, tauP)
%NLAYERRECTANGULARWEIGHTSANDNODES Summary of this function goes here
%   Todo: Write documentation
%
%   Designation for each dimension:
%       1: Integration variable tau
%       2: Mode matrix columns
%       3: Mode matrix rows
%       4: Subscript i of integration parameters I_i(n,m,p,q)
%       5: Integration variable psi (temporarily)

arguments
    O
    tauP(:, 1)
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
am(1, 1, :, 1) = O.modes(:, 1) * pi ./ O.a;
bn(1, 1, :, 1) = O.modes(:, 2) * pi ./ O.b;
ap(1, :, 1, 1) = O.modes(:, 1) * pi ./ O.a;
bq(1, :, 1, 1) = O.modes(:, 2) * pi ./ O.b;

%% Helper Functions
f = @(xi, eta) (O.a.^2 .* O.b.^2 / 16) ...
    .* sinc((0.5/pi)*O.a * (xi - am)) ...
    .* sinc((0.5/pi)*O.a * (xi - ap)) ...
    .* sinc((0.5/pi)*O.b * (eta - bn)) ...
    .* sinc((0.5/pi)*O.b * (eta - bq)) ...
    ./ ((am + xi) .* (ap + xi) ...
    .* (bn + eta) .* (bq + eta));

% f = @(xi, eta) (O.a.^2 .* O.b.^2 / 16) .* (0.25./am./ap) ...
%     .* (sinc((0.5/pi)*O.a * (xi - am)) + sinc((0.5/pi)*O.a * (xi + am))) ...
%     .* (sinc((0.5/pi)*O.a * (xi - ap)) + sinc((0.5/pi)*O.a * (xi + ap))) ...
%     .* sinc((0.5/pi)*O.b * (eta - bn)) .* sinc((0.5/pi)*O.b * (eta - bq)) ...
%     ./ ((bn + eta) .* (bq + eta));

%% Compute Integrals Over Psi at All Values of Tau
[psi(1, 1, 1, 1, :), weightsPsi(1, 1, 1, 1, :)] = ...
    O.gaussLegendre(O.integralPointsPsi, 0, 0.5*pi);
xi = tau .* cos(psi);
eta = tau .* sin(psi);

fsc = 4*sum(weightsPsi .* tau.^3 .* sin(psi).^2 .* cos(psi).^2 .* f(xi, eta), 5);
fss = 4*sum(weightsPsi .* tau.^3 .* sin(psi).^4 .* f(xi, eta), 5);
fcc = 4*sum(weightsPsi .* tau.^3 .* cos(psi).^4 .* f(xi, eta), 5);

%% Compute Integrand Values
integralE = zeros(size(fsc)) + zeros(1, 1, 1, 4);
integralH = zeros(size(integralE));

integralE(:, :, :, 1) = weights .* (bn .* ap .* (4/pi.^2) .* fsc);
integralE(:, :, :, 2) = weights .* (am .* ap .* (4/pi.^2) .* fss);
integralE(:, :, :, 3) = weights .* (bn .* bq .* (4/pi.^2) .* fcc);
integralE(:, :, :, 4) = weights .* (am .* bq .* (4/pi.^2) .* fsc);

integralH(:, :, :, 1) = -weights .* (bn .* ap .* (4/pi.^2) .* fsc);
integralH(:, :, :, 2) =  weights .* (am .* ap .* (4/pi.^2) .* fsc);
integralH(:, :, :, 3) =  weights .* (bn .* bq .* (4/pi.^2) .* fsc);
integralH(:, :, :, 4) = -weights .* (am .* bq .* (4/pi.^2) .* fsc);

%% Fix Nans Caused by Singularities
integralE(tauP == 0, :, :, :) = 0;
integralH(tauP == 0, :, :, :) = 0;

integralE(tauP == 1, :, :, :) = 0;
integralH(tauP == 1, :, :, :) = 0;

end
