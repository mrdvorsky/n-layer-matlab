function [Ex, Ey] = getModeFields(self, x, y)
%Return the fields defined by this "waveguideMode" object.
% This function returns the tangential fields over the waveguide aperture.
%
% Author: Matt Dvorsky

arguments
    self waveguideMode;

    x(:, 1);
    y(1, :);
end

%% Calculate kx, ky, kr, and kphi
[kx, ky] = fftCoordinates(x, y, ApplyFftShift=true);
kr = hypot(kx, ky);
kphi = atan2(ky, kx);

%% Calculate Fields
ExHat = self.ExSpec(kx, ky, kr, kphi);
EyHat = self.EySpec(kx, ky, kr, kphi);

ExHat(:, 1) = 0;
EyHat(:, 1) = 0;
ExHat(1, :) = 0;
EyHat(1, :) = 0;

intScaleFactor = numel(kr) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
    ./ (2*pi);

Ex = fftshift(ifft2(ifftshift(ExHat))) * intScaleFactor;
Ey = fftshift(ifft2(ifftshift(EyHat))) * intScaleFactor;

end

