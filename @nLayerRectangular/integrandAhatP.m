function [Ahat_kRhoP] = integrandAhatP(O, kRhoP, k0, er, ur, thk)
%INTEGRANDA1 Integrand for the matrix A. Integrate over [0, 1].
% This function computes the integrand used to compute the integral of
% A. The matrix A can thus be computed using the following example.
%
% Example Usage:
%   A = integral(@(kRhoP) O.integrandAhatP(kRhoP, k0, er, ur, thk), ...
%       0, 1, "ArrayValued", "true");
%
% After integrating this function over [0, 1], the output can be used 
% along with the outputs of the constructFrequencyMultipliers(...) function
% and the constructMatrixEquation(...) function to calculate S11 for a 
% rectangular waveguide. See the documentation of "computeA1" for more
% details.
%
% This function is meant to be used inside of an adaptive integration
% method, and shouldn't be used in fixed point techniques. This is because
% of the optimization used and explained in the first section below.
%
% This function is used in the "computeA" function. See documentation for
% more details.
%
% Inputs:
%   kRhoP: Column vector of coordinates at which to evaluate AhatP(kRhoP).
%   k0 - Free-space wavenumber (must have compatible dimensions with kRhoP).
%   er - Array of complex relative permittivities for each layer (must have
%       compatible dimensions with k0).
%   ur - Array of complex relative permeabilities for each layer (must have
%       compatible dimensions with k0).
%   thk - Array of thicknesses for each layer (must have compatible 
%       dimensions with er and ur). The last element of thk should have a 
%       value of inf in the infinite halfspace case.
% Output:
%   Ahat_kRhoP - Ah(kRhoP) calculated for each kRhoP and k0. The dimensions
%       of Ahat_kRhoP will be numel(kRhoP) by O.numModes by O.numModes by ...,
%       where the last dimensions are based on k0.
%
% Author: Matt Dvorsky

%% First Pass Optimization
% The initial pass of the adaptive integration routine used in
% "computeA" uses the same nodes every time. We can thus take advantage
% of precomputed values for A1_EH instead of performing the same
% interpolation every time. We can tell if this is the first pass by
% comparing the size of kRhoP to the size of the precomputed A1_EH. Note
% that this size is asserted to be an odd integer, and that only the first
% pass of an adaptive integration routine can have an odd value for the
% size of kRhoP. Thus, there is no need to check the values of kRhoP.
%
% The values of O.init_kRho, O.init_A1_E, and O.init_A1_H are computed
% in the "recomputeInterpolants" member function.
if numel(kRhoP) == numel(O.init_kRho)
    [GammaH, GammaE] = O.computeGamma0(O.init_kRho, k0, er, ur, thk);
    Ahat_kRhoP = GammaH .* O.init_A1_H + GammaE .* O.init_A1_E;
    return;
end

%% General Case (Linear Interpolation)
kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
[GammaH, GammaE] = O.computeGamma0(kRho, k0, er, ur, thk);

% Get indices and mixing factors for linear interpolation
fracInd = kRhoP * (O.interpolationPoints_kRho - 1) + 1;
intInd = floor(fracInd);
m = fracInd - intInd;

% Perform linear interpolation
vLower = O.A1_EH(intInd, :, :, :);
vHigher = O.A1_EH(intInd + 1, :, :, :);
interp_A1_EH = vLower + m .* (vHigher - vLower);

Ahat_kRhoP = GammaE.*interp_A1_EH(:, :, :, 1) ...
    + GammaH.*interp_A1_EH(:, :, :, 2);

end


