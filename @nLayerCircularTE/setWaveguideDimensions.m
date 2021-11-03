function [] = setWaveguideDimensions(O, r)
%SETWAVEGUIDEDIMENSIONS Set circular waveguide radius.
% Calling this functions sets radius of the circular waveguide to the
% specified value. The default unit is mm, however, the units must be
% changed if not using the default value of the speed of light ("c"), which
% is defined in the nLayerForward class.
%
% After calling this function, the "recomputeInterpolants" function should
% be called before calling "calculate".
%
% Example Usage:
%   NL = nLayerCircularTE(...);
%   NL.setWaveguideDimensions(r);
%   NL.recomputeInterpolants();
%
% Inputs:
%   r - New value of O.r (waveguide radius), in mm.
%
% Author: Matt Dvorsky

arguments
    O;
    r(1, 1) {mustBeNumeric, mustBePositive};
end

O.r = r;

end

