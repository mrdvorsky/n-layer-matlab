function [varargout] = solveStructure(O, NL, f, gam, options)
%SOLVESTRUCTURE Perform curve fitting to solve for structure parameters.
% This function takes triplets of nLayerForward objects, frequency vectors,
% and measurements, and tries to find the missing structure parameters of
% er, ur, thk using curve fitting to minimize the rms values of
% 'NL.calculate(f, er, ur, thk) - gam' for each triplet.
% 
% Each NL object can be different, as long as 'NL.calculate(f, ...)'
% returns an array with the same size as 'gam'. The first dimension of
% 'gam' must have a size match the number of elements in 'f'.
%
% Example Usage (simple case, 2 layer infinite halfspace):
%   NLsolver.setInitialValues(Er=[1, 4-1j], Thk=[10, inf]);
%   NLsolver.setLayersToSolve(Erp=[2], Erpp=[2]);
%   [Params, Gamma, Uncert] = NLsolver.solveStructure(NL, f, gam);
%
% Example Usage (for multi-band open-ended measurements):
%   NLsolver.setInitialValues(Er=[1, 4-1j], Thk=[1, 10]);
%   NLsolver.setLayersToSolve(Erp=[2], Erpp=[2], Thk=[1, 2]);
%   [Params, Gamma, Uncert] = NLsolver.solveStructure(...
%       NL1, f1, gam1, ...
%       NL2, f2, gam2);
%
% Inputs:
%   NL (Repeating) - A valid nLayerForward object.
%   f (Repeating) - Vector of frequencies to pass to NL.
%   gam (Repeating) - Measurements to fit. Size must match the output size
%       of 'NL.calculate(f, er, ur, thk)'.
%
% Outputs:
%   Parameters - Struct containing the structure parameters (er, ur, thk),
%       and simulated measurements (gam) for each input set.
%   Gamma - Cell array of simulated measurements (i.e., the value of
%       'NL.calculate(f, er, ur, thk)').
%   Uncertainty - Struct containing the calculated output parmeter
%       uncertainties for each input set.
%
% Named Arguments:
%   NoiseStdMin (0.001) - Minimum uncertainty value to assume for the
%       measurement data when calculating the Uncertainty struct. If the
%       RMS difference between the fit and the measurements is less than
%       NoiseStdMin, NoiseStdMin will be used instead.
%
% Author: Matt Dvorsky

arguments
    O;
end

arguments(Repeating)
    NL(1, 1) {mustBeA(NL, "nLayerForward")};
    f(:, 1) {mustBeNonempty};
    gam {mustBeCorrectGamSize(f, gam)};
end

arguments
    options.NoiseStdMin(1, 1) {mustBeNonnegative} = 0.01;
end

%% Perform Curve Fitting Using 'solveStructureMultiple'
inputParams = [repmat({O}, 1, numel(NL)); NL; f; gam];
[varargout{1:nargout}] = nLayerInverse.solveStructureMultiple(...
    inputParams{:}, NoiseStdMin=options.NoiseStdMin);

if nargout >= 1
    varargout{1} = varargout{1}{1};     % Get first element of Params.
end
if nargout >= 3
    varargout{3} = varargout{3}{1};     % Get first element of Uncert.
end

end


function mustBeCorrectGamSize(f, gam)
    if iscell(f)    % Fix MATLAB bug.
        currentInd = find(cellfun(@(x) numel(x) > 0, f), 1, "last");
        f = f{currentInd};
    end
    if numel(f) ~= size(gam, 1)
        throwAsCaller(MException("nLayerInverse:mustBeCorrectGamSize", ...
            "First dimension of the measurements array must have size " + ...
            "equal to the number of frequencies (%d).", numel(f)));
    end
end
