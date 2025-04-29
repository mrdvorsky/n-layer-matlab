clc;
clear;
close all;

%% Inputs
% NL = nLayerCircular(0, 5, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
% f = linspace(32, 40, 801);

NL = nLayerRectangular(3, 2, waveguideBand="X");
f = linspace(8.2, 12.4, 801);

numPointsX = 1001;
numPointsY = 1001;

z = 3*linspace(20, -50, 200);

%% Calculate x, y, kx, and ky
sizeX = 1.2 * max([NL.modeStructs.ApertureWidth]);
sizeY = sizeX;
x(:, 1) = sizeX * linspace(-0.5, 0.5, numPointsX);
y(1, :) = sizeY * linspace(-0.5, 0.5, numPointsY);

[kx, ky] = fftCoordinates(x, y);

%% Calculate Fields
for ii = flip(1:numel(NL.modeStructs))
    ExHat = NL.modeStructs(ii).ExSpec(kx, ky, hypot(kx, ky), atan2(ky, kx));
    EyHat = NL.modeStructs(ii).EySpec(kx, ky, hypot(kx, ky), atan2(ky, kx));

    scaleFactor = numel(ExHat) * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
        ./ (2*pi);
    Ex(:, :, ii) = ifftshift(ifft2(ExHat)) * scaleFactor + zeros(numel(kx), numel(ky));
    Ey(:, :, ii) = ifftshift(ifft2(EyHat)) * scaleFactor + zeros(numel(kx), numel(ky));
end

%% Plot Fields
figure;
ModeImg = showImage(x, y, abs(Ex(:, :, 1)));
clim([0, inf]);

S = NL.calculate(f, 4 - 0.004j, [], 2);

%% Plot Layer Fields
[ExSpec2, EySpec2, kx2, ky2, xL] = get_spectrums(NL);
[Ex2, Ey2] = get_fields(z, f(end), {1}, {1}, {10}, ...
    zeros(NL.numModes, 1), ExSpec2, EySpec2, ...
    kx2, ky2, [NL.modeStructs.CutoffWavenumber]);

figure;
Img = showImage(xL, z, Ey2, DisplayFormat="Imag");
colormap colormapplusminusabs;

%% Viewer
NL2 = copy(NL);
NL2.receiveModeIndices = 1:NL.numModes;

figure;
nLayerViewer({4 - 0.004j}, {}, {2}, NL, f, ...
    UpdateFunction=@(f, er, ur, thk) updateFun(...
    f, er, ur, thk, ...
    NL2, ...
    ModeImg, Ex, Ey, ...
    Img, z, ExSpec2, EySpec2, kx2, ky2, NL.mode_kc0));



%% Helper
function updateFun(f, er, ur, thk, NL, ModeImg, Ex, Ey, Img, z, ExSpec2, EySpec2, kx2, ky2, beta_c)
    S_single(1, 1, :) = NL.calculate(f, er, ur, thk);
    exc(1, 1, :) = 1:numel(S_single) == 1;
    Ex = sum(Ex .* (S_single + exc), 3);
    Ey = sum(Ey .* (S_single + exc), 3);
    phaseMax = 0.5 * angle(mean((Ex(:) + Ey(:)).^2));
    u = real(Ex .* exp(-1j .* phaseMax));
    v = real(Ey .* exp(-1j .* phaseMax));
    ModeImg.CData = hypot(u, v).';
    title(ModeImg.Parent, sprintf("f = %g GHz, Phase = %g deg", f, rad2deg(phaseMax)));

    % [~, Ey2] = get_fields(z, f, er, ur, thk, ...
    %     S_single, ExSpec2, EySpec2, ...
    %     kx2, ky2, beta_c);
    % Img.CData = real(Ey2 .* exp(-1j .* phaseMax)).';
end

function [ExSpec, EySpec, kx, ky, x] = get_spectrums(NL)
    numPointsX = 1024;
    numPointsY = 1024;

    sizeX = 10 * max([NL.modeStructs.ApertureWidth]);
    sizeY = sizeX;
    x(:, 1) = sizeX * linspace(-0.5, 0.5, numPointsX);
    y(1, :) = sizeY * linspace(-0.5, 0.5, numPointsY);
    
    [kx, ky] = fftCoordinates(x, y);

    for ii = flip(1:numel(NL.modeStructs))
        scaleFactor = numel(kx) * numel(ky) ...
            * abs(kx(1) - kx(2)) * abs(ky(1) - ky(2)) ...
            ./ (2*pi);
        ExSpec(:, :, ii) = NL.modeStructs(ii).ExSpec(kx, ky, ...
            hypot(kx, ky), atan2(ky, kx)) .* scaleFactor ...
            + zeros(numel(kx), numel(ky));;
        EySpec(:, :, ii) = NL.modeStructs(ii).EySpec(kx, ky, ...
            hypot(kx, ky), atan2(ky, kx)) .* scaleFactor ...
            + zeros(numel(kx), numel(ky));;
    end
end

function [Ex, Ey] = get_fields(z, f, er, ur, thk, gam, ExSpecAll, EySpecAll, kx, ky, beta_c)
    gam = reshape(gam, 1, 1, []);
    gam(1) = gam(1);
    z = reshape(z, 1, 1, 1, []);

    k0 = 2*pi*f/299.792;
    kz = sqrt(k0.^2.*er{1}.*ur{1} - kx.^2 - ky.^2);
    beta(1, 1, :) = sqrt(k0.^2 - beta_c.^2);
    kz = complex(real(kz), -abs(imag(kz)));
    beta = complex(real(beta), -abs(imag(beta)));

    ExSpec = sum(ExSpecAll.*gam, 3) + ExSpecAll(:, :, 1);
    EySpec = sum(EySpecAll.*gam, 3) + EySpecAll(:, :, 1);
    Ex2 = fftshift(ifft(sum(ExSpec .* exp(1j*(z(z<0)).*kz), 2), [], 1), 1);
    Ey2 = fftshift(ifft(sum(EySpec .* exp(1j*(z(z<0)).*kz), 2), [], 1), 1);

    Ex1 = sum(gam .* fftshift(ifft(sum(ExSpecAll, 2), [], 1), 1) ...
        .* exp(-1j*(beta).*z(z>=0)), 3);
    Ey1 = sum(gam .* fftshift(ifft(sum(EySpecAll, 2), [], 1), 1) ...
        .* exp(-1j*(beta).*z(z>=0)), 3);
    Ex1 = Ex1 + sum(fftshift(ifft(sum(ExSpecAll(:, :, 1), 2), [], 1), 1) ...
        .* exp(1j*beta(1).*z(z>=0)), 3);
    Ey1 = Ey1 + sum(fftshift(ifft(sum(EySpecAll(:, :, 1), 2), [], 1), 1) ...
        .* exp(1j*beta(1).*z(z>=0)), 3);

    Ex = cat(2, Ex1(:, :), Ex2(:, :));
    Ey = cat(2, Ey1(:, :), Ey2(:, :));
end



