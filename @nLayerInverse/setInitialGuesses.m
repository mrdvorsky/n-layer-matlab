function [] = setInitialGuesses(O, options)
%SETINITIALGUESS Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
    options.ErGuess(1, :) {mustBePositive, mustBeFinite, ...
        mustBeGreaterThanOrEqual(options.ErGuess, 1)};
    options.ErpGuess(1, :) {mustBePositive, mustBeFinite};
    options.ThkGuess(1, :) {mustBePositive, mustBeFinite};
end

%% Assign Initial Guesses
O.erGuess =  options.ErGuess;
O.erpGuess = options.ErpGuess;
O.thkGuess = options.ThkGuess;

end

