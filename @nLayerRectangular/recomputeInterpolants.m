function [] = recomputeInterpolants(O)
%RECOMPUTEINTERPOLANTS Recompute interpolation functions, structs, etc.
% This function should be called anytime any critical parameter is changed,
% and must be called before calling calculate. This function is
% automatically called after creating an nLayerRectangular object.
%
% List of critical parameters:
%   a;
%   b;
%   speedOfLight; (Inherited from nLayerForward)
%   modesTE;
%   interpolationPointsTau;
%   integralPointsTauFixed;
%   integralInitialSegmentCount;
%   integralPointsPsi;
%
% Example Usage:
%   NL = nLayerRectangular(...);
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

%% Update TM Modes
% Every valid TM mode having a corresponding TE mode will be considered.
O.modesTM = O.modesTE(find(O.modesTE(:, 2) > 0), :);
O.numModes = size(O.modesTE, 1) + size(O.modesTM, 1);

% Scale factor for change of variables between tau and tauP
O.integralScaleFactor = pi*pi ./ O.a;

%% Compute A1 and b1 at various values of tauP
% Compute integrands for I_i^(e) and I_i^(h) at each value of tauP.
tauP(:, 1) = linspace(0, 1, O.interpolationPointsTau);
[integrandE, integrandH] = O.computeIntegrandEH(tauP);

% Construct matrix equation from integrands. Note that a permutation is
% first performed, as "constructMatrixEquation" expects the first 3
% dimensions to be m, n, and i.
% Note that A2 does not depend on I_i^(e) or I_i^(h), and thus is only
% computed once.
[A1_E, A2] = O.constructMatrixEquation(permute(integrandE, [2, 3, 4, 1]));
[A1_H, ~] = O.constructMatrixEquation(permute(integrandH, [2, 3, 4, 1]));

% Combine A1_E and A1_H into one array for faster interpolation.
% Also, undo the previous permutation.
A1_EH = ipermute(cat(3, A1_E, A1_H), [2, 3, 4, 1]);

% Store computed matrices
O.A1_EH = A1_EH;
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
[integrandE, integrandH] = O.computeIntegrandEH(tauP);
[A1_E, ~] = O.constructMatrixEquation(permute(integrandE, [2, 3, 4, 1]));
[A1_H, ~] = O.constructMatrixEquation(permute(integrandH, [2, 3, 4, 1]));

A1_E = ipermute(A1_E, [2, 3, 4, 1]);
A1_H = ipermute(A1_H, [2, 3, 4, 1]);

% Store computed matrices. Also, precompute tau using tauP.
O.fixed_tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
O.fixed_A1_E = A1_E .* weights;
O.fixed_A1_H = A1_H .* weights;
O.fixed_errA1_E = A1_E .* errWeights;
O.fixed_errA1_H = A1_H .* errWeights;

%% First Integration Pass Precomputation
% The initial pass of the adaptive integral algorithm always uses the same
% nodes (i.e., the same evaluation coordinates of tau). Thus, we can 
% precompute the values of the integrand for A1 at those coordinates.
% These are used in the "integrandA1" function.
[tauP, ~, ~] = O.gaussKronrod(...
    O.integralInitialSegmentCount, 0, 1);

% The procedure here is the same as in the previous section.
[integrandE, integrandH] = O.computeIntegrandEH(tauP);
[A1_E, ~] = O.constructMatrixEquation(permute(integrandE, [2, 3, 4, 1]));
[A1_H, ~] = O.constructMatrixEquation(permute(integrandH, [2, 3, 4, 1]));

A1_E = ipermute(A1_E, [2, 3, 4, 1]);
A1_H = ipermute(A1_H, [2, 3, 4, 1]);

% Store computed matrices. Also, precompute tau using tauP.
O.init_tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
O.init_A1_E = A1_E;
O.init_A1_H = A1_H;

end

