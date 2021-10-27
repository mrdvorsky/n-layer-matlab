function [] = setWaveguideDimensions(O, a, b)
%SETWAVEGUIDEDIMENSIONS Set waveguide broad and narrow dimensions.
% Calling this functions sets the broad and narrow dimensions, O.a and O.b,
% of the rectangular waveguide, to the specified values. The default unit
% is mm, however, the units must be changed if not using the default value
% of the speed of light ("c"), which is defined in the nLayerForward class.
%
% After calling this function, the "recomputeInterpolants" function should
% be called before calling "calculate".
%
% Example Usage:
%   NL = nLayerRectangular(...);
%   NL.setWaveguideDimensions(a, b);
%   NL.recomputeInterpolants();
%
% Inputs:
%   a - New value of O.a (waveguide broad dimension), in mm.
%   b - New value of O.b (waveguide narrow dimension), in mm.
%
% Author: Matt Dvorsky

arguments
    O;
    a(1, 1) {mustBeNumeric, mustBePositive};
    b(1, 1) {mustBeNumeric, mustBePositive};
end

O.a = a;
O.b = b;

end

