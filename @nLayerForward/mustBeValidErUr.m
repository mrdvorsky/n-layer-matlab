function [] = mustBeValidErUr(er)
%MUSTBEVALIDERUR Throws an error if input is not a valid er or ur.
% Throws an error if er is not a valid relative permittivity or
% permeability. More specifically, it will throw an error if the real part
% is less than 1 or if the imaginary part is positive.
%
% Example Usage:
%   arguments
%       er(1, :) {nLayerForward.mustBeValidErUr};
%   end
%
% Author: Matt Dvorsky

arguments
    er;
end

%% Check Argument
mustBeFinite(er);
if any(real(er) < 1, "all")
    throwAsCaller(MException("nLayer:mustBeValidErUr", ...
        "The real part must be 1 or more."));
end
if any(imag(er) > 0, "all")
    throwAsCaller(MException("nLayer:mustBeValidErUr", ...
        "The imaginary part must be nonnegative."));
end

end

