function [erOut, urOut, thkOut] = extractStructure(O, x, f)
%EXTRACTSTRUCTURE Converts the linearized structure (x) to er, ur, and thk.
% This is a helper function that takes the linearized structure parameter
% guess (x) passed by a curve fitting solver and returns the full parameter
% structures er, ur, and thk. The linearized structure (x) must have been
% created using the same nLayerInverse settings.
%
% Author: Matt Dvorsky

arguments
    O;
    x(:, 1);
    f(:, 1);
end

%% Create Default Structure
erp  =  real(O.initialValue_er);
erpp = -imag(O.initialValue_er);
urp  =  real(O.initialValue_ur);
urpp = -imag(O.initialValue_ur);
thk  =  O.initialValue_thk;

%% Fill in Structure
xInds = cumsum([0, length(O.layersToSolve_erp), length(O.layersToSolve_erpp), ...
    length(O.layersToSolve_urp), length(O.layersToSolve_urpp), ...
    length(O.layersToSolve_thk)]);

x_erp  = x(xInds(1) + 1:xInds(2));
x_erpp = x(xInds(2) + 1:xInds(3));
x_urp  = x(xInds(3) + 1:xInds(4));
x_urpp = x(xInds(4) + 1:xInds(5));
x_thk  = x(xInds(5) + 1:xInds(6));

erp(1,  O.layersToSolve_erp)  = x_erp.';
erpp(1, O.layersToSolve_erpp) = x_erpp.';
urp(1,  O.layersToSolve_urp)  = x_urp.';
urpp(1, O.layersToSolve_urpp) = x_urpp.';
thk(1,  O.layersToSolve_thk)  = x_thk.';

%% Construct Outputs
erOut = complex(erp, -erpp);
urOut = complex(urp, -urpp);
thkOut = thk;

end

