function [] = setWaveguideDimensions(O, a, b)
%SETWAVEGUIDEDIMENSIONS Set waveguide a and b dimensions.
% Calling this functions sets O.a and O.b to the specified values. The 
% default unit is mm, however, the units must be changed if not using the
% default value of the speed of light ("c"), which is defined in the
% nLayerForward class.
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
%   a - New value of O.a
%   b - New value of O.b
%
% Author: Matt Dvorsky

O.a = a;
O.b = b;

end

