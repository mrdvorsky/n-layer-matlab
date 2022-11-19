function [gamError, gamErrorComplex] = calculateError(NLsolver, x, NL, f, gamActual, options)
%CALCULATEERROR Helper function to calculate error in a given curve fit iteration.
% This function takes the linearized structure parameter guess (x)
% specified by the curve fitting function and calculates the error vector
% (or mse) between the measurements (gam) and the simulated measurements
% from the forward solvers (NL). If the 'VectorOutput' named argument is
% set to false, will return mse of the error vector.
%
% Author: Matt Dvorsky

arguments
    NLsolver;
    x(:, 1);
    NL(:, 1) cell;
    f(:, 1) cell;
    gamActual(:, 1) cell;
    options.VectorOutput = true;
end

%% Construct Multilayer Structure
[er, ur, thk] = NLsolver.extractStructure(x, f);

%% Calculate Gamma
gam = cell(length(NL), 1);
try
    for ii = 1:length(NL)
        gam{ii} = NL{ii}.calculate(f{ii}, er, ur, thk);
    end
catch ex
    error("Failed to evaluate structure because: %s\n%s", ex.message, ...
        NLsolver.printStructureParameters(er, ur, thk, Title="Failed to Converge"));
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
    gamError = sum(gamError.^2);
end

end

