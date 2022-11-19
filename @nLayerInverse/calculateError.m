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
    NLsolver(:, 1) cell;
    x(:, 1);
    NL(:, 1) cell;
    f(:, 1) cell;
    gamActual(:, 1) cell;
    options.VectorOutput = true;
end

%% Calculate Gamma for Each Structure
gam = cell(numel(NLsolver), 1);
try
    for ii = 1:numel(NLsolver)
        [er, ur, thk] = NLsolver{ii}.extractStructure(x, f);
        gam{ii} = NL{ii}.calculate(f{ii}, er, ur, thk);
    end
catch ex
    error("Failed to evaluate structure because: %s\n%s", ex.message, ...
        NLsolver{ii}.printStructureParameters(er, ur, thk, Title="Failed to Converge"));
end

%% Calculate Error
gamErrorComplex = cell(numel(NLsolver), 1);
for ii = 1:numel(NLsolver)
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

