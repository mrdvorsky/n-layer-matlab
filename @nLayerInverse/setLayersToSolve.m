function [] = setLayersToSolve(O, options)
%SETLAYERSTOSOLVE Sets the layer indices to solve for each variable.
% This function is simply a setter for the layersToSolve_{...} parameters
% that does some bounds checking. Any unspecified parameter indices will
% be left unchanged.
%
% Example Usage:
%   NLsolver.setLayersToSolve(Erp=[1, 2], Erpp=[2]);
%   NLsolver.setLayersToSolve(Urp=[1], Urpp=[2], Thk=[1, 2]);
%
% Author: Matt Dvorsky

arguments
    O;
    options.Erp (:, 1) {mustBeValidLayerInd(options.Erp,  O)} = O.layersToSolve_erp;
    options.Erpp(:, 1) {mustBeValidLayerInd(options.Erpp, O)} = O.layersToSolve_erpp;
    options.Urp (:, 1) {mustBeValidLayerInd(options.Urp,  O)} = O.layersToSolve_urp;
    options.Urpp(:, 1) {mustBeValidLayerInd(options.Urpp, O)} = O.layersToSolve_urpp;
    options.Thk (:, 1) {mustBeValidLayerInd(options.Thk,  O)} = O.layersToSolve_thk;
end

%% Assign Layer Indices
O.layersToSolve_erp  = options.Erp;
O.layersToSolve_erpp = options.Erpp;
O.layersToSolve_urp  = options.Urp;
O.layersToSolve_urpp = options.Urpp;
O.layersToSolve_thk  = options.Thk;

end

function mustBeValidLayerInd(inds, O)
    if any(inds > O.layerCount) || numel(inds) ~= numel(unique(inds)) ...
            || any(inds <= 0) || any(inds ~= floor(inds))
        throwAsCaller(MException("nLayerInverse:mustBeValidLayerInd", ...
            "The parameters 'layersToSolve_{er, ur, thk}' must " + ...
            "consist of unique positive integers no greater than " + ...
            "the layer count (%d).", O.layerCount));
    end
end
