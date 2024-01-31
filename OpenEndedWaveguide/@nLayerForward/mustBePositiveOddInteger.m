function [] = mustBePositiveOddInteger(num)
%MUSTBEPOSITIVEODDINTEGER Throws an error when not a positive odd integer.
%
% Example Usage:
%   arguments
%       posOddInt {nLayerForward.mustBePositiveOddInteger};
%   end
%
% Author: Matt Dvorsky

arguments
    num;
end

%% Check Argument
mustBeInteger(num);
mustBePositive(num);
if mod(num, 2) ~= 1
    throwAsCaller(MException("MATLAB:mustBePositiveOddInteger", ...
        "Value must be a positive odd integer."));
end

end

