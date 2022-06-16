function [AHat] = integrandAhat(O, kRhoP, k0, er, ur, thk)
%INTEGRANDAHAT Integrand for the matrix A. Integrate over [0, 1].
% This function computes the integrand used to compute the integral of
% the matrix A. While this function can used to compute A directly using
% the example below, it is recommended to use the "computeA" function
% instead, which uses this internally. See documentation of "computeA" for
% more details.
%
% Example Usage:
%   A = integral(@(kRhoP) O.integrandAhatP(kRhoP, k0, er, ur, thk), ...
%       0, 1, "ArrayValued", "true");
%
% This function is meant to be used inside of an adaptive integration
% method, and shouldn't be used in fixed point techniques. This is because
% of the optimization used and explained in the first section below.
%
% Inputs:
%   kRhoP: Column vector of coordinates at which to evaluate Ahat(kRhoP).
%   k0 - Free-space wavenumber (must have compatible dimensions with kRhoP).
%   er - Array of complex relative permittivities for each layer (must have
%       compatible dimensions with k0).
%   ur - Array of complex relative permeabilities for each layer (must have
%       compatible dimensions with k0).
%   thk - Array of thicknesses for each layer (must have compatible 
%       dimensions with er and ur). The last element of thk should have a 
%       value of inf in the infinite halfspace case.
% Output:
%   Ahat - Ahat(kRhoP) calculated for each kRhoP and k0. The dimensions
%       of Ahat_kRhoP will be numel(kRhoP) by O.numModes by O.numModes by ...,
%       where the last dimensions are based on k0.
%
% Author: Matt Dvorsky

%% First Pass Optimization
% The initial pass of the adaptive integration routine used in
% "computeA" uses the same nodes every time. We can thus take advantage
% of precomputed values for AhHat and AeHat instead of performing the same
% interpolation every time. We can tell if this is the first pass by
% comparing the size of kRhoP to the size of the precomputed AhHat. Note
% that this size is asserted to be an odd integer, and that only the first
% pass of an adaptive integration routine can have an odd value for the
% size of kRhoP. Thus, there is no need to check that the values of kRhoP
% match the precomputed values.
%
% The values of O.init_kRho, O.init_AhHat, and O.init_AeHat are computed
% in the "recomputeInterpolants" member function.
if numel(kRhoP) == numel(O.init_kRho)
    [Gamma0h, Gamma0e] = O.computeGamma0(O.init_kRho, k0, er, ur, thk);
    AHat = Gamma0h .* O.init_AhHat + Gamma0e .* O.init_AeHat;
    return;
end

%% General Case (Linear Interpolation)
kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
[Gamma0h, Gamma0e] = O.computeGamma0(kRho, k0, er, ur, thk);

% Get indices and mixing factors for linear interpolation
fracInd = kRhoP * (O.interpolationPoints_kRho - 1) + 1;
intInd = floor(fracInd);
mixingFactor = fracInd - intInd;

% Perform linear interpolation
vLower = O.table_AheHat(intInd, :, :, :);
vHigher = O.table_AheHat(intInd + 1, :, :, :);
interpolated_AheHat = vLower + mixingFactor .* (vHigher - vLower);

AHat = Gamma0h .* interpolated_AheHat(:, :, :, 1) ...
    + Gamma0e .* interpolated_AheHat(:, :, :, 2);

end


