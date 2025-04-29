function [] = mustBeErUrCallable(erFun)
%Throws an error if input is not a valid er or ur object.
% Throws an error if "erFun" is not a valid relative permittivity or
% permeability object. A valid object will either be a numeric scalar or a
% function that accepts a frequency argument and returns a value.
%
% Example Usage:
%   arguments
%       erFun(1, 1) {nLayer.mustBeErUrCallable};
%       erScalar(1, 1) {nLayer.mustBeErUrCallable};
%   end
%
% Author: Matt Dvorsky

arguments
    erFun(:, 1);
end

%% Check Input
% Numeric input.
if isnumeric(erFun)
    return;
end

% Single function handle.
if ~iscell(erFun)
    try
        mustBeCallable(erFun, {1}, "freq");
    catch ex
        throwAsCaller(MException("nLayer:mustBeErUrCallable", ...
            sprintf("Value must be one of the following:" + ...
            "\n\t- A scalar value." + ...
            "\n\t- A callable object with call syntax 'erur = erurFun(f)'." + ...
            "\n\t- An array of values." + ...
            "\n\t- A cell array, where each element is a callable or scalar.")));
    end
    return;
end

% Cell array of function handles or scalar values.
for ii = 1:numel(erFun)
    if isnumeric(erFun{ii}) && numel(erFun{ii}) == 1
        continue;
    end

    try
        mustBeCallable(erFun{ii}, {1}, "freq");
    catch ex
        throwAsCaller(MException("nLayer:mustBeErUrCallable", ...
            sprintf("Value must be one of the following:" + ...
            "\n\t- A scalar value." + ...
            "\n\t- A callable object with call syntax 'erur = erurFun(f)'." + ...
            "\n\t- An array of values." + ...
            "\n\t- A cell array, where each element is a callable or scalar.")));
    end
end

end

