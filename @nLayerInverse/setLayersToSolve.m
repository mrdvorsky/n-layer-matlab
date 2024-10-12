function [] = setLayersToSolve(O, options)
%SETLAYERSTOSOLVE Sets the layer indices to solve for each variable.
% This function is simply a setter for the layersToSolve_{...} parameters
% that does some bounds checking. Any unspecified parameter indices will
% be left unchanged.
%
% Example Usage:
%   NLsolver.setLayersToSolve(Er=[1], Thk=[1]);
%   NLsolver.setLayersToSolve(Ur=[1, 2], Thk=[1, 2]);
%   NLsolver.setLayersToSolve(Erp=[1], Erpp=[1, 2]);    % Separate erp and erpp.
%
% Author: Matt Dvorsky

arguments
    O;
    options.Er  (1, :) = [];
    options.Ur  (1, :) = [];
    options.Thk (1, :) = O.layersToSolve_thk;

    options.Erp (1, :) = O.layersToSolve_erp;
    options.Erpp(1, :) = O.layersToSolve_erpp;
    options.Urp (1, :) = O.layersToSolve_urp;
    options.Urpp(1, :) = O.layersToSolve_urpp;
end

%% Assign Layer Indices
if ~isempty(options.Er)
    options.Erp = options.Er;
    options.Erpp = options.Er;
end
if ~isempty(options.Ur)
    options.Urp = options.Ur;
    options.Urpp = options.Ur;
end

O.layersToSolve_erp  = options.Erp;
O.layersToSolve_erpp = options.Erpp;
O.layersToSolve_urp  = options.Urp;
O.layersToSolve_urpp = options.Urpp;
O.layersToSolve_thk  = options.Thk;

end

