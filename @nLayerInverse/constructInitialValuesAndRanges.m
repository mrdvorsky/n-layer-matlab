function [xInitial, xMin, xMax] = constructInitialValuesAndRanges(O)
%CONSTRUCTINITIALVALUESANDRANGES Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
end

%% Check Value Bounds
if ~isempty(O.erLayersToSolve)
    if any(O.erInitialValue(O.erLayersToSolve) < O.erRange(1, :) ...
            | O.erInitialValue(O.erLayersToSolve) > O.erRange(2, :))
        error("An element of 'erInitialValue' is outside the range specified by 'erRange'.");
    end
end

if ~isempty(O.erpLayersToSolve)
    if any(O.erpInitialValue(O.erpLayersToSolve) < O.erpRange(1, :) ...
            | O.erpInitialValue(O.erpLayersToSolve) > O.erpRange(2, :))
        error("An element of 'erpInitialValue' is outside the range specified by 'erpRange'.");
    end
end

if ~isempty(O.thkLayersToSolve)
    if any(O.thkInitialValue(O.thkLayersToSolve) < O.thkRange(1, :) ...
            | O.thkInitialValue(O.thkLayersToSolve) > O.thkRange(2, :))
        error("An element of 'thkInitialValue' is outside the range specified by 'thkRange'.");
    end
end

if any(~isfinite(O.thkInitialValue(O.thkLayersToSolve)))
    error("Last layer thickness cannot be infinite if it is being solved.");
end

%% Construct Initial Values and Ranges
erInitialValue  = O.erInitialValue( 1, O.erLayersToSolve);
erpInitialValue = O.erpInitialValue(1, O.erpLayersToSolve);
thkInitialValue = O.thkInitialValue(1, O.thkLayersToSolve);

erMin  = (O.erRange( 1, O.erLayersToSolve));
erpMin = (O.erpRange(1, O.erpLayersToSolve));
thkMin = O.thkRange(1, O.thkLayersToSolve);

erMax  = (O.erRange( 2, O.erLayersToSolve));
erpMax = (O.erpRange(2, O.erpLayersToSolve));
thkMax = O.thkRange(2, O.thkLayersToSolve);

%% Assemble Output
% Transformation of er and erp is done to improve convergence. If this is
% changed, make sure to make the corresponding change in the
% "extractStructure" function and the min and max values above.
xInitial = [(erInitialValue), (erpInitialValue), thkInitialValue].';
xMin   = [erMin,   erpMin,   thkMin].';
xMax   = [erMax,   erpMax,   thkMax].';

end

