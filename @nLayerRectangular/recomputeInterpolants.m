function [] = recomputeInterpolants(O)
%NLAYERRECTCOMPUTEINTEGRALSTRUCT Summary of this function goes here
%   Detailed explanation goes here

%% Compute A1 and b1 at various values of tauP
% Scale factor for change of variables between tau and tauP
O.integralScaleFactor = pi*pi ./ O.a;

% Compute integrands for I_i^e and I_i^h at each value of tauP
tauP(:, 1) = linspace(0, 1, O.interpolationPointsTau);
[integralE, integralH] = O.computeIntegralEH(tauP);

% Construct matrix equation from integrands. Note that a permutation is
% first performed, as constructMatrixEquation() expects the first 3
% dimensions to be m, n, and i.
[A1_E, A2, b1_E, b2] = O.constructMatrixEquation(permute(integralE, [2, 3, 4, 1]));
[A1_H, ~, b1_H, ~] = O.constructMatrixEquation(permute(integralH, [2, 3, 4, 1]));

% Combine A1_E, b1_E, A1_H, b1_H into one array for faster interpolation.
% Also, undo the previous permutation.
A1b1_EH = ipermute(cat(3, [A1_E, b1_E], [A1_H, b1_H]), [2, 3, 4, 1]);

% Store computed matrices
O.A1b1_EH = A1b1_EH;
O.A2 = A2;
O.b2 = b2;

%% Lossy Structure Weights and Nodes
% For lossy structures, generally no adaptive meshing is needed. In those
% cases we can use a precomputed set of weights and nodes, instead of
% computing them on the fly. This is generally 3 to 4 times faster than
% when using the adaptive integration.
[tauP, weights, errWeights] = O.gaussKronrod(...
    round(O.integralPointsTauLossy / 15), 0, 1);

[integralE, integralH] = O.computeIntegralEH(tauP);
[A1_E, ~, b1_E, ~] = O.constructMatrixEquation(permute(integralE, [2, 3, 4, 1]));
[A1_H, ~, b1_H, ~] = O.constructMatrixEquation(permute(integralH, [2, 3, 4, 1]));

A1b1_E = ipermute([A1_E, b1_E], [2, 3, 4, 1]);
A1b1_H = ipermute([A1_H, b1_H], [2, 3, 4, 1]);

tau = O.integralScaleFactor * (1 - tauP) ./ tauP;

% Store computed matrices
O.lossy_tau = tau;
O.lossy_A1b1_E = A1b1_E .* weights;
O.lossy_A1b1_H = A1b1_H .* weights;
O.lossy_errA1b1_E = A1b1_E .* errWeights;
O.lossy_errA1b1_H = A1b1_H .* errWeights;

O.init_A1b1_E = A1b1_E;
O.init_A1b1_H = A1b1_H;

end

