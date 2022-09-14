function [erOut, urOut, thkOut] = extractStructure(O, x, f)
%EXTRACTSTRUCTURE Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
    x(:, 1);
    f(:, 1);
end

%% Create Default Structure
er = O.erInitialValue;
erp = O.erpInitialValue;
ur = O.urInitialValue;
urp = O.urpInitialValue;
thk = O.thkInitialValue;

%% Fill in Structure
xInds = cumsum([0, length(O.erLayersToSolve), length(O.erpLayersToSolve), ...
    length(O.urLayersToSolve), length(O.urpLayersToSolve), length(O.thkLayersToSolve)]);

x_er  = x(xInds(1) + 1:xInds(2));
x_erp = x(xInds(2) + 1:xInds(3));
x_ur = x(xInds(3) + 1:xInds(4));
x_urp = x(xInds(4) + 1:xInds(5));
x_thk = x(xInds(5) + 1:xInds(6));

% Transformation of er and erp is done to improve convergence. If this is
% changed, make sure to make the corresponding change in the
% "constructInitialValuesAndRanges" function.
er(1, O.erLayersToSolve) = (x_er).';
erp(1, O.erpLayersToSolve) = (x_erp).';
ur(1, O.urLayersToSolve) = (x_ur).';
urp(1, O.urpLayersToSolve) = (x_urp).';
thk(1, O.thkLayersToSolve) = x_thk.';

%% Construct Outputs
erOut = complex(er, -erp);
urOut = complex(ur, -urp);
thkOut = thk;

end

