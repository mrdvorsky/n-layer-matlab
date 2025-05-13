function [] = validateModeAmplitude(self, options)
%Check that the fields have a complex-amplitude of 1.
%
% Author: Matt Dvorsky

arguments
    self waveguideMode;

    options.NumIntegralPointsRho(1, 1) {mustBeInteger, mustBePositive} = 8000;
    options.NumIntegralPointsPhi(1, 1) {mustBeInteger, mustBePositive} = 200;

    options.AbsTol(1, 1) {mustBePositive} = 1e-6;
end

%% Calculate Mode Self Cross Product
modeAmp = waveguideMode.modeCrossProduct(self, self, ...
    NumIntegralPointsRho=options.NumIntegralPointsRho, ...
    NumIntegralPointsPhi=options.NumIntegralPointsPhi);

%% Check Complex-Amplitude Value
if abs(modeAmp - 1) > options.AbsTol

    error("waveguideMode:incorrectAmplitide", ...
        "The amplitude of mode of mode '%s' is (%.16g + j%.16g), " + ...
        "but it should be (1 - j0), within a tolerance of (%g).", ...
        self.modeLabel, ...
        real(modeAmp), imag(modeAmp), ...
        options.AbsTol);
end

end




