function [outputArg1,outputArg2] = mustBeErUrCallable(erFun)
%MUSTBECALLABLEERUR Throws an error if input is not a valid er or ur object.
% Throws an error if "erFun" is not a valid relative permittivity or
% permeability object. A valid object will either be a numeric scalar or a
% function that accepts a frequency argument and returns a value.
%
% Example Usage:
%   arguments
%       erFun(1, 1) {nLayer.mustBeValidErUr};
%       erScaler(1, 1) {nLayer.mustBeValidErUr};
%   end
%
% Author: Matt Dvorsky

arguments
    erFun(1, 1);
end

%% Check Input
if isnumeric(erFun)
    return;
end

mustBeCallable(erFun, {1}, "freq");

end

