function [gamError, gamErrorComplex] = calculateError(O, x, NL, f, gamActual, options)
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
try
    gam = NL.calculate(f, er, ur, thk);
catch ex
    error("Failed to evaluate structure because: %s\n%s", ex.message, ...
        O.printStructureParameters(er, ur, thk, Title="Falied to Converge"));
end

%% Calculate Error
gamErrorComplex = gam(:) - gamActual(:);
gamError = [real(gamErrorComplex); imag(gamErrorComplex)] ...
    ./ sqrt(numel(gamErrorComplex));

if ~options.VectorOutput
    gamError = (sum(gamError.^2));
end

end

