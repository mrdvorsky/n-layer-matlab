function [gamError, gamErrorComplex] = calculateError(O, x, NL, f, gamActual, options)
%CALCULATEERROR Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
    x(:, 1);
    NL(:, 1) cell;
    f(:, 1) cell;
    gamActual(:, 1) cell;
    options.VectorOutput = true;
end

%% Construct Multilayer Structure
[er, ur, thk] = O.extractStructure(x, f);

%% Calculate Gamma
gam = cell(length(NL), 1);
try
    for ii = 1:length(NL)
        gam{ii} = NL{ii}.calculate(f{ii}, er, ur, thk);
    end
catch ex
    error("Failed to evaluate structure because: %s\n%s", ex.message, ...
        O.printStructureParameters(er, ur, thk, Title="Falied to Converge"));
end

%% Calculate Error
gamErrorComplex = cell(length(NL), 1);
for ii = 1:length(NL)
    gamErrorComplex{ii} = gam{ii}(:) - gamActual{ii}(:);
end

%% Unify Output
gamErrorOutput = cat(1, gamErrorComplex{:});
gamError = [real(gamErrorOutput); imag(gamErrorOutput)] ...
    ./ sqrt(numel(gamErrorOutput));

if ~options.VectorOutput
    gamError = (sum(gamError.^2));
end

end

