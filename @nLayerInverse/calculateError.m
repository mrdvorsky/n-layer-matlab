function [gamError] = calculateError(O, x, NL, f, gamActual, options)
%CALCULATEERROR Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
    x(:, 1);
    NL;
    f(:, 1);
    gamActual;
    options.VectorOutput = true;
end

%% Construct Multilayer Structure
[er, ur, thk] = O.extractStructure(x, f);

%% Calculate Gamma
gam = NL.calculate(f, er, ur, thk);

%% Calculate Error
gamErrorComplex = gam(:) - gamActual(:);
gamError = [real(gamErrorComplex); imag(gamErrorComplex)];

if ~options.VectorOutput
    gamError = rms(gamError);
end

end

