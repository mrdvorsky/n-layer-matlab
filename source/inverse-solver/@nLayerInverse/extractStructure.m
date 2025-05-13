function [erOut, urOut, thkOut] = extractStructure(self, x, f)
%Converts the linearized structure (x) to er, ur, and thk.
% This is a helper function that takes the linearized structure parameter
% guess (x) passed by a curve fitting solver and returns the full parameter
% structures er, ur, and thk. The linearized structure (x) must have been
% created using the same nLayerInverse settings.
%
% Author: Matt Dvorsky

arguments
    self nLayerInverse;

    x(:, 1);
    f(:, 1);
end

%% Create Default Structure
erp  =  real(self.initialValue_er);
erpp = -imag(self.initialValue_er);
urp  =  real(self.initialValue_ur);
urpp = -imag(self.initialValue_ur);
thk  =  self.initialValue_thk;

%% Fill in Structure
xInds = cumsum([0, length(self.layersToSolve_erp), length(self.layersToSolve_erpp), ...
    length(self.layersToSolve_urp), length(self.layersToSolve_urpp), ...
    length(self.layersToSolve_thk)]);

x_erp  = x(xInds(1) + 1:xInds(2));
x_erpp = x(xInds(2) + 1:xInds(3));
x_urp  = x(xInds(3) + 1:xInds(4));
x_urpp = x(xInds(4) + 1:xInds(5));
x_thk  = x(xInds(5) + 1:xInds(6));

erp(1,  self.layersToSolve_erp)  = x_erp.';
erpp(1, self.layersToSolve_erpp) = x_erpp.';
urp(1,  self.layersToSolve_urp)  = x_urp.';
urpp(1, self.layersToSolve_urpp) = x_urpp.';
thk(1,  self.layersToSolve_thk)  = x_thk.';

%% Construct Outputs
erOut = complex(erp, -erpp);
urOut = complex(urp, -urpp);
thkOut = thk;

end

