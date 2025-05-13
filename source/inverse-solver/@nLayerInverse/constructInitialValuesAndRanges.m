function [xInitial, xMin, xMax, xA, xb, xAeq, xbeq] = constructInitialValuesAndRanges(self)
%Linearizes initial values and ranges for optimizer.
% Performs the conversion from structure parameters (er, ur, thk) to the
% linearized real vector parameters taken in by the curve fitting
% functions. This process makes use of the layersToSolve_{...} parameters
% to determine which layers to be solved.
%
% This function also checks that the initial values being solved are
% within the ranges specified by the rangeMin_{...} and rangeMax_{...}
% parameters.
%
% Author: Matt Dvorsky

arguments
    self nLayerInverse;
end

%% Construct Initial Values and Ranges
erp  =  real(self.initialValue_er( self.layersToSolve_erp));
erpp = -imag(self.initialValue_er( self.layersToSolve_erpp));
urp  =  real(self.initialValue_ur( self.layersToSolve_urp));
urpp = -imag(self.initialValue_ur( self.layersToSolve_urpp));
thk  =       self.initialValue_thk(self.layersToSolve_thk);

erpMin  = self.rangeMin_erp( self.layersToSolve_erp);
erppMin = self.rangeMin_erpp(self.layersToSolve_erpp);
urpMin  = self.rangeMin_urp( self.layersToSolve_urp);
urppMin = self.rangeMin_urpp(self.layersToSolve_urpp);
thkMin  = self.rangeMin_thk( self.layersToSolve_thk);

erpMax  = self.rangeMax_erp( self.layersToSolve_erp);
erppMax = self.rangeMax_erpp(self.layersToSolve_erpp);
urpMax  = self.rangeMax_urp( self.layersToSolve_urp);
urppMax = self.rangeMax_urpp(self.layersToSolve_urpp);
thkMax  = self.rangeMax_thk( self.layersToSolve_thk);

%% Apply Thickness Constraints
thk_Aeq = self.constraints_thk_Aeq;
thk_beq = thk_Aeq * self.initialValue_thk(:);

thk_A = self.constraints_thk_A;
thk_b = self.constraints_thk_b;

%% Assemble Output
xInitial = [erp, erpp, urp, urpp, thk].';
xMin = [erpMin, erppMin, urpMin, urppMin, thkMin].';
xMax = [erpMax, erppMax, urpMax, urppMax, thkMax].';

% Linear thickness constraints (xAeq*x == xbeq and xA*x <= xb).
xAeq = zeros(numel(thk_beq), numel(xMax));
xA = zeros(numel(thk_b), numel(xMax));
xbeq = thk_beq;
xb = thk_b;

xAeq(:, end-numel(thkMax)+1:end) = thk_Aeq(:, self.layersToSolve_thk);
xA(:,   end-numel(thkMax)+1:end) = thk_A(:,   self.layersToSolve_thk);

%% Check Value Bounds
if any(xInitial < xMin, "all") || any(xInitial > xMax, "all")
    error("One or more structure parameters (er, ur, thk) that " + ...
        "are being solved for is outside the range specified.");
end

if any(~isfinite(thk))
    error("Last layer thickness cannot be infinite if it is being solved.");
end

end


