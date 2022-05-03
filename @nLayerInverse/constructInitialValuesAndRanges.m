function [xInitial, xMin, xMax] = constructInitialValuesAndRanges(O)
%CONSTRUCTINITIALVALUESANDRANGES Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
end

%% Check Guess Bounds
if any(O.erGuess < O.erRange(1, :) | O.erGuess > O.erRange(2, :))
    error("An element of 'erGuess' is outside the range specified by 'erRange'.");
end

if any(O.erpGuess < O.erpRange(1, :) | O.erpGuess > O.erpRange(2, :))
    error("An element of 'erpGuess' is outside the range specified by 'erpRange'.");
end

if any(O.thkGuess < O.thkRange(1, :) | O.thkGuess > O.thkRange(2, :))
    error("An element of 'thkGuess' is outside the range specified by 'thkRange'.");
end

if any(~isfinite(O.thkGuess(O.thkLayersToSolve)))
    error("Last layer thickness cannot be infinite if it is being solved.");
end

%% Construct Guesses and Ranges
erGuess  = O.erGuess( 1, O.erLayersToSolve);
erpGuess = O.erpGuess(1, O.erpLayersToSolve);
thkGuess = O.thkGuess(1, O.thkLayersToSolve);

erMin  = O.erRange( 1, O.erLayersToSolve);
erpMin = O.erpRange(1, O.erpLayersToSolve);
thkMin = O.thkRange(1, O.thkLayersToSolve);

erMax  = O.erRange( 2, O.erLayersToSolve);
erpMax = O.erpRange(2, O.erpLayersToSolve);
thkMax = O.thkRange(2, O.thkLayersToSolve);

%% Assemble Output
xInitial = [erGuess, erpGuess, thkGuess].';
xMin   = [erMin,   erpMin,   thkMin].';
xMax   = [erMax,   erpMax,   thkMax].';

end

