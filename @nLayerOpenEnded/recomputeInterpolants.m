function [] = recomputeInterpolants(O, options)
%RECOMPUTEINTERPOLANTS Recompute interpolation functions, structs, etc.
% This function should be called anytime any critical parameter is changed,
% and must be called before calling calculate. This function is
% automatically called after creating an nLayerRectangular object.
%
% List of critical parameters:
%   waveguideA;
%   waveguideB;
%   speedOfLight;
%   modesTE;
%   interpolationPoints_kRho;
%   integralPointsFixed_kRho;
%   integralPoints_kPhi;
%   integralInitialSegmentCount;
%
% Example Usage:
%   NL = nLayerRectangular(...);
%   NL.*criticalParameter* = *newValue*;
%   NL.recomputeInterpolants();
%
% Author: Matt Dvorsky

arguments
    O;

    options.SpecEx_TE(:, 1) {mustBeA(options.SpecEx_TE, "cell")} = {};
    options.SpecEy_TE(:, 1) {mustBeA(options.SpecEy_TE, "cell")} = {};
    options.SpecEx_TM(:, 1) {mustBeA(options.SpecEx_TM, "cell")} = {};
    options.SpecEy_TM(:, 1) {mustBeA(options.SpecEy_TM, "cell")} = {};
    options.SpecEx_Hybrid(:, 1) {mustBeA(options.SpecEx_Hybrid, "cell")} = {};
    options.SpecEy_Hybrid(:, 1) {mustBeA(options.SpecEy_Hybrid, "cell")} = {};

    options.ModeSymmetryX {mustBeMember(options.ModeSymmetryX, ...
        ["None", "Even", "Odd"])} = "None";
    options.ModeSymmetryY {mustBeMember(options.ModeSymmetryY, ...
        ["None", "Even", "Odd"])} = "None";
    options.ModeSymmetryAxial {mustBeMember(options.ModeSymmetryAxial, ...
        ["None", "TE", "TM"])} = "None";

    options.IntegralScaleFactor(1, 1) {mustBePositive} = 1;
end

%% Check Input Mode Sizes
[options.SpecEx_TE, options.SpecEy_TE] = ...
    checkSpectrumSizes(options.SpecEx_TE, options.SpecEy_TE);
[options.SpecEx_TM, options.SpecEy_TM] = ...
    checkSpectrumSizes(options.SpecEx_TM, options.SpecEy_TM);
[options.SpecEx_Hybrid, options.SpecEy_Hybrid] = ...
    checkSpectrumSizes(options.SpecEx_Hybrid, options.SpecEy_Hybrid);

%% Update Mode Counts
O.numModes_TE = numel(options.SpecEx_TE);
O.numModes_TM = numel(options.SpecEx_TM);
O.numModes_Hybrid = numel(options.SpecEx_Hybrid);
O.numModes = O.numModes_TE + O.numModes_TM + O.numModes_Hybrid;

%% Update Integral Scale Factor
O.integralScaleFactor = options.IntegralScaleFactor;

%% Compute Ahat at all possible values of kRhoP
kRhoP(:, 1) = linspace(0, 1, O.interpolationPoints_kRho);

% Compute AhHat and AeHat at kRhoP coordinates.
[AhHat, AeHat] = O.computeAhat(kRhoP, options);

% Combine AeHat and AhHat into one array for faster interpolation.
% Store in a table that is used in the "integrandAhat" function.
O.table_AheHat = cat(4, AhHat, AeHat);

%% Fixed Point Integration Weights and Nodes
% For lossy structures, generally no adaptive meshing is needed. In those
% cases we can use a precomputed set of weights and nodes, instead of
% computing them on the fly. This is generally 3 to 4 times faster than
% when using the adaptive integration.
[kRhoP, weights, errWeights] = O.fejer2(O.integralPointsFixed_kRho, 0, 1);

% Compute AhHat and AeHat at kRhoP coordinates.
[AhHat, AeHat] = O.computeAhat(kRhoP, options);

% Store computed matrices. Also, precompute kRho using kRhoP. These are
% used in the "computeA" function.
O.fixed_kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
O.fixed_AhHat = AhHat .* weights;
O.fixed_AeHat = AeHat .* weights;
O.fixed_errorAhHat = AhHat .* errWeights;
O.fixed_errorAeHat = AeHat .* errWeights;

%% First Integration Pass Precomputation
% The initial pass of the adaptive integral algorithm always uses the same
% nodes (i.e., the same evaluation coordinates of kRho). Thus, we can
% precompute the values of the integrand for A at those coordinates,
% instead of needing to perform an interpolation. These are used in the
% "integrandAhat" function.
[kRhoP, ~, ~] = O.gaussKronrod(...
    O.integralInitialSegmentCount, 0, 1);

% Compute AhHat and AeHat at kRhoP coordinates.
[AhHat, AeHat] = O.computeAhat(kRhoP, options);

% Store computed matrices. Also, precompute kRho using kRhoP. These are
% used in the "integrandAhat" function.
O.init_kRho = O.integralScaleFactor * (1 - kRhoP) ./ kRhoP;
O.init_AhHat = AhHat;
O.init_AeHat = AeHat;

end



%% Input Argument Checking
function [specEx, specEy] = checkSpectrumSizes(specEx, specEy)

if isempty(specEx) && ~isempty(specEy)
    specEx = repmat({@(~, ~, ~, ~) 0}, numel(specEy), 1);
end
if isempty(specEy) && ~isempty(specEx)
    specEy = repmat({@(~, ~, ~, ~) 0}, numel(specEx), 1);
end
if numel(specEx) ~= numel(specEy)
    error("Mode spectrum arguments for Ex and Ey must " + ...
        "be cell arrays of the same size or one must be empty.");
end

end

