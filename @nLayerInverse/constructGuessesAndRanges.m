function [xGuess, xMin, xMax] = constructGuessesAndRanges(O)
%CONSTRUCTGUESSESANDRANGES Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
end

%% Check Guess Bounds
if any(O.erGuess < O.erRange(1, :) | O.erGuess > O.erRange(2, :))
    error("erGuess is outside the range specified by erRange.");
end

if any(O.erpGuess < O.erpRange(1, :) | O.erpGuess > O.erpRange(2, :))
    error("erpGuess is outside the range specified by erpRange.");
end

if any(O.thkGuess < O.thkRange(1, :) | O.thkGuess > O.thkRange(2, :))
    error("thkGuess is outside the range specified by thkRange.");
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
xGuess = [erGuess, erpGuess, thkGuess].';
xMin   = [erMin,   erpMin,   thkMin].';
xMax   = [erMax,   erpMax,   thkMax].';

end

