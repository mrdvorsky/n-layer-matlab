function [] = setWaveguideDimensions(O, waveguideA, waveguideB)
%SETWAVEGUIDEDIMENSIONS Set waveguide broad and narrow dimensions.
% Calling this functions sets the broad and narrow dimensions, O.a and O.b,
% of the rectangular waveguide, to the specified values. The default unit
% is mm, however, the units must be changed if not using the default value
% of SpeedOfLight, which is defined in the nLayerForward class.
%
% After calling this function, the "recomputeInterpolants" function should
% be called before calling "calculate".
%
% Example Usage:
%   NL = nLayerRectangular(...);
%   NL.setWaveguideDimensions(waveguideA, waveguideB);
%   NL.recomputeInterpolants();
%
% Inputs:
%   waveguideA - New value of O.waveguideA (waveguide broad dimension).
%   waveguideB - New value of O.waveguideB (waveguide narrow dimension).
%
% Author: Matt Dvorsky

arguments
    O;
    waveguideA(1, 1) {mustBeNumeric, mustBePositive};
    waveguideB(1, 1) {mustBeNumeric, mustBePositive};
end

%% Set Dimensions
O.waveguideA = waveguideA;
O.waveguideB = waveguideB;

end

