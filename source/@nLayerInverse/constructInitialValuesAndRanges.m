function [xInitial, xMin, xMax, xA, xb, xAeq, xbeq] = constructInitialValuesAndRanges(O)
%CONSTRUCTINITIALVALUESANDRANGES Linearizes initial values and ranges for optimizer.
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
    O;
end

%% Construct Initial Values and Ranges
erp  =  real(O.initialValue_er( O.layersToSolve_erp));
erpp = -imag(O.initialValue_er( O.layersToSolve_erpp));
urp  =  real(O.initialValue_ur( O.layersToSolve_urp));
urpp = -imag(O.initialValue_ur( O.layersToSolve_urpp));
thk  =       O.initialValue_thk(O.layersToSolve_thk);

erpMin  = O.rangeMin_erp( O.layersToSolve_erp);
erppMin = O.rangeMin_erpp(O.layersToSolve_erpp);
urpMin  = O.rangeMin_urp( O.layersToSolve_urp);
urppMin = O.rangeMin_urpp(O.layersToSolve_urpp);
thkMin  = O.rangeMin_thk( O.layersToSolve_thk);

erpMax  = O.rangeMax_erp( O.layersToSolve_erp);
erppMax = O.rangeMax_erpp(O.layersToSolve_erpp);
urpMax  = O.rangeMax_urp( O.layersToSolve_urp);
urppMax = O.rangeMax_urpp(O.layersToSolve_urpp);
thkMax  = O.rangeMax_thk( O.layersToSolve_thk);

%% Apply Thickness Constraints
thk_Aeq = O.constraints_thk_Aeq;
thk_beq = thk_Aeq * O.initialValue_thk(:);

thk_A = O.constraints_thk_A;
thk_b = O.constraints_thk_b;

%% Assemble Output
xInitial = [erp, erpp, urp, urpp, thk].';
xMin = [erpMin, erppMin, urpMin, urppMin, thkMin].';
xMax = [erpMax, erppMax, urpMax, urppMax, thkMax].';

% Linear thickness constraints (xAeq*x == xbeq and xA*x <= xb).
xAeq = zeros(numel(thk_beq), numel(xMax));
xA = zeros(numel(thk_b), numel(xMax));
xbeq = thk_beq;
xb = thk_b;

xAeq(:, end-numel(thkMax)+1:end) = thk_Aeq(:, O.layersToSolve_thk);
xA(:,   end-numel(thkMax)+1:end) = thk_A(:,   O.layersToSolve_thk);

%% Check Value Bounds
if any(xInitial < xMin, "all") || any(xInitial > xMax, "all")
    error("One or more structure parameters (er, ur, thk) that " + ...
        "are being solved for is outside the range specified.");
end

if any(~isfinite(thk))
    error("Last layer thickness cannot be infinite if it is being solved.");
end

end


