function [] = recomputeInterpolants(O)
%RECOMPUTEINTERPOLANTS Recompute interpolation functions, structs, etc.
% This function should be called anytime any critical parameter is changed,
% and must be called before calling calculate. This function is
% automatically called after creating an nLayerCircularTE object.
%
% List of critical parameters:
%   r;
%   speedOfLight; (Inherited from nLayerForward)
%   modesTE;
%   interpolationPointsTau;
%   integralPointsTauFixed;
%   integralInitialSegmentCount;
%   integralPointsPsi;
%
% Example Usage:
%   NL = nLayerCircularTE(...);
%   NL.*criticalParameter* = *newValue*;
%   NL.recomputeInterpolants();
%
% Author: Matt Dvorsky

%% Check For Parameter Validity
% This initial segment count for adaptive integration should be odd so
% that the "integrandA1" function can use that fact to determine when the
% first pass occurs. See "integrandA1" for more details.
if mod(O.integralInitialSegmentCount, 2) == 0 % Check if even
    error(strcat("Parameter 'integralInitialSegmentCount' must be an ", ...
        "odd integer (current value: %d)."), O.integralInitialSegmentCount);
end

%% Calculate Mode Cutoffs
O.numModes = numel(O.modesTE);

O.modeCutoffs = zeros(O.numModes, 1);
for ii = 1:O.numModes
    % Use zeros of 1st order bessel function to calculate mode cutoff
    % wavenumbers. A reasonable approximation for the nth zero is 
    % 0.785 + pi*n. This is used as the initial guess.
    O.modeCutoffs(ii) = fzero(@(x) besselj(1, x), ...
        0.785 + pi*O.modesTE(ii)) ./ O.r;
end

% Scale factor for change of variables between tau and tauP
O.integralScaleFactor = (2*pi) ./ O.r;

%% Compute A1 and b1 at various values of tauP
% Compute integrands for I_i^(e) and I_i^(h) at each value of tauP.
tauP(:, 1) = linspace(0, 1, O.interpolationPointsTau);
integrandH = O.computeIntegrandH(tauP);

% Construct matrix equation from integrands. Note that a permutation is
% first performed, as "constructMatrixEquation" expects the first 3
% dimensions to be m, n, and i.
% Note that A2 and b2 do not depend on I^(h), and thus are only computed
% once. Also note that b1 is ignored since it is always equal to the
% negative of the first column of A1.
[A1_H, A2] = O.constructMatrixEquation(permute(integrandH, [2, 3, 4, 1]));

% Undo the previous permutation.
A1_H = ipermute(A1_H, [2, 3, 4, 1]);

% Store computed matrices
O.A1_H = A1_H;
O.A2 = A2;

%% Fixed Point Integration Weights and Nodes
% For lossy structures, generally no adaptive meshing is needed. In those
% cases we can use a precomputed set of weights and nodes, instead of
% computing them on the fly. This is generally 3 to 4 times faster than
% when using the adaptive integration.
[tauP, weights, errWeights] = O.fejer2(O.integralPointsTauFixed, 0, 1);

% The procedure here is almost exactly the same as in the previous section,
% except there is no need to recompute A2 and b2, and A1_E and A2_E are
% kept separate.
integrandH = O.computeIntegrandH(tauP);
[A1_H, ~] = O.constructMatrixEquation(permute(integrandH, [2, 3, 4, 1]));

A1_H = ipermute(A1_H, [2, 3, 4, 1]);

% Store computed matrices. Also, precompute tau using tauP.
O.fixed_tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
O.fixed_A1_H = A1_H .* weights;
O.fixed_errA1_H = A1_H .* errWeights;

%% First Integration Pass Precomputation
% The initial pass of the adaptive integral algorithm always uses the same
% nodes (i.e., the same evaluation coordinates of tau). Thus, we can 
% precompute the values of the integrand for A1 at those coordinates.
% These are used in the "integrandA1" function..
[tauP, ~, ~] = O.gaussKronrod(...
    O.integralInitialSegmentCount, 0, 1);

% The procedure here is the same as in the previous section.
integrandH = O.computeIntegrandH(tauP);
[A1_H, ~] = O.constructMatrixEquation(permute(integrandH, [2, 3, 4, 1]));

A1_H = ipermute(A1_H, [2, 3, 4, 1]);

% Store computed matrices. Also, precompute tau using tauP.
O.init_tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
O.init_A1_H = A1_H;

end

