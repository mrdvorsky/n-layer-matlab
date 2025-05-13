function [] = validateModeSymmetry(self, options)
%Check that the fields have the symmetry specified.
%
% Author: Matt Dvorsky

arguments
    self waveguideMode;

    options.NumSamplesRho(1, 1) {mustBeInteger, mustBePositive} = 100;
    options.NumSamplesPhi(1, 1) {mustBeInteger, mustBePositive} = 100;
    options.AbsTol(1, 1) {mustBePositive} = 1e-6;
end

%% Sample Coordinates
Nrho = options.NumSamplesRho;
Nphi = options.NumSamplesPhi;

krMax = 10 * 2*pi ./ self.apertureSize;

kr(:, 1)  = fejer2(Nrho, 0, krMax);
kphi(1, :, 1) = fejer2(Nphi, 0, 0.5*pi);
kphi(1, :, 2) = fejer2(Nphi, pi, 0.5*pi);
kphi(1, :, 3) = fejer2(Nphi, pi, 1.5*pi);
kphi(1, :, 4) = fejer2(Nphi, 0, -0.5*pi);

%% Sample Fields
kx = kr .* cos(kphi);
ky = kr .* sin(kphi);

Wh = self.WhSpec(kx, ky, kr, kphi) + 0*kx;
We = self.WeSpec(kx, ky, kr, kphi) + 0*kx;

%% Check Symmetry
xPEC = max(abs(Wh(:, :, [1, 4]) + conj(Wh(:, :, [2, 3]))), [], "all") ...
    + max(abs(We(:, :, [1, 4]) - conj(We(:, :, [2, 3]))), [], "all");
yPEC = max(abs(Wh(:, :, [1, 2]) + conj(Wh(:, :, [4, 3]))), [], "all") ...
    + max(abs(We(:, :, [1, 2]) - conj(We(:, :, [4, 3]))), [], "all");
xPMC = max(abs(Wh(:, :, [1, 4]) - conj(Wh(:, :, [2, 3]))), [], "all") ...
    + max(abs(We(:, :, [1, 4]) + conj(We(:, :, [2, 3]))), [], "all");
yPMC = max(abs(Wh(:, :, [1, 2]) - conj(Wh(:, :, [4, 3]))), [], "all") ...
    + max(abs(We(:, :, [1, 2]) + conj(We(:, :, [4, 3]))), [], "all");

symmetryX = "None";
if xPEC < options.AbsTol
    symmetryX = "PEC";
elseif xPMC < options.AbsTol
    symmetryX = "PMC";
end

symmetryY = "None";
if yPEC < options.AbsTol
    symmetryY = "PEC";
elseif yPMC < options.AbsTol
    symmetryY = "PMC";
end

%% Axial Symmetry
axialTE_std = max(std(Wh(:, :), 1, 2));
axialTM_std = max(std(We(:, :), 1, 2));
axialTE_max = max(abs(Wh(:)));
axialTM_max = max(abs(We(:)));

symmetryAxial = "None";
if (axialTE_std < options.AbsTol) && (axialTM_max < options.AbsTol)
    symmetryAxial = "TE";
    symmetryX = "PEC";
    symmetryY = "PEC";
elseif (axialTM_std < options.AbsTol) && (axialTE_max < options.AbsTol)
    symmetryAxial = "TM";
    symmetryX = "PMC";
    symmetryY = "PMC";
end

%% Compare
if (symmetryX ~= self.symmetryX) ...
        || (symmetryY ~= self.symmetryY) ...
        || (symmetryAxial ~= self.symmetryAxial)
    error("waveguideMode:symmetryMismatch", ...
        "Symmetry of mode of mode '%s' is specified as " + ...
        "(x='%s', y='%s', axial='%s'), but is actually " + ...
        "(x='%s', y='%s', axial='%s').", ...
        self.modeLabel, ...
        self.symmetryX, self.symmetryY, self.symmetryAxial, ...
        symmetryX, symmetryY, symmetryAxial);
end

end




