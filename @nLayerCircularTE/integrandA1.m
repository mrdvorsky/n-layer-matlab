function [A1_tauP] = integrandA1(O, tauP, k0, er, ur, thk)
%INTEGRANDA1 Integrand for the matrix A1. Integrate over [0, 1].
% This function computes the integrand used to compute the integral of
% A1. The matrix A1 can thus be computed using the following example.
%
% Example Usage:
%   A1 = integral(@(tauP) O.integrandA1(tauP, k0, er, ur, thk), 0, 1, ...
%       "ArrayValued", "true");
%
% After integrating this function over [0, 1], the output can be used 
% along with the outputs of the constructFrequencyMultipliers(...) function
% and the constructMatrixEquation(...) function to calculate S11 for a 
% circular TE0n mode waveguide. See the documentation of "computeA1" for
% more details.
%
% This function is meant to be used inside of an adaptive integration
% method, and shouldn't be used in fixed point techniques. This is because
% of the optimization used and explained in the first section below.
%
% This function is used in the "computeA1" function. See documentation for
% more details.
%
% Inputs:
%   tauP: Column vector of coordinates at which to evaluate A1(tauP).
%   k0 - Free-space wavenumber (must have compatible dimensions with tau).
%   er - Array of complex relative permittivities for each layer (must have
%       compatible dimensions with k0).
%   ur - Array of complex relative permeabilities for each layer (must have
%       compatible dimensions with k0).
%   thk - Array of thicknesses for each layer (must have compatible 
%       dimensions with er and ur). The last element of thk should have a 
%       value of inf in the infinite halfspace case.
% Output:
%   A1_tauP - A1(tauP) calculated for each tauP and k0. The dimensions of
%       A1_tauP will be numel(tauP) by O.numModes by O.numModes by ...,
%       where the last dimensions are based on k0.
%
% Author: Matt Dvorsky

%% First Pass Optimization
% The initial pass of the adaptive integration routine used in
% "calculateA1" uses the same nodes every time. We can thus take advantage
% of precomputed values for A1_H instead of performing the same
% interpolation every time. We can tell if this is the first pass by
% comparing the size of tauP to the size of the precomputed A1_H. Note
% that this size is asserted to be an odd integer, and that only the first
% pass of an adaptive integration routine can have an odd value for the
% size of tauP. Thus, there is no need to check the values of tauP.
%
% The values of O.init_tau and O.init_A1_H are computed in the
% "recomputeInterpolants" member function.
if numel(tauP) == numel(O.init_tau)
    specH = O.multilayerSpectrumH(O.init_tau, k0, er, ur, thk);
    A1_tauP = specH .* O.init_A1_H;
    return;
end

%% General Case (Linear Interpolation)
tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
specH = O.multilayerSpectrumH(tau, k0, er, ur, thk);

% Get indices and mixing factors for linear interpolation
fracInd = tauP * (O.interpolationPointsTau - 1) + 1;
intInd = floor(fracInd);
m = fracInd - intInd;

% Perform linear interpolation
vLower = O.A1_H(intInd, :, :, :);
vHigher = O.A1_H(intInd + 1, :, :, :);
interp_A1_H = vLower + m .* (vHigher - vLower);

A1_tauP = specH .* interp_A1_H;

end


