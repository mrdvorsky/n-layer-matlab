clc;
clear;
close all;

%% Inputs
NL = nLayerRectangular(1, 2, waveguideBand="Ka");

ExFun = NL.modeStructs(end).ExSpec;
EyFun = NL.modeStructs(end).EySpec;

%% Coordinates
[kx, ky] = fftCoordinates(10*linspace(-5, 5, 1000), 10*linspace(-2.5, 2.5, 1000), ApplyFftShift=true);

kphi = atan2(ky, kx);
kr = hypot(kx, ky);

scale = (kx(2) - kx(1)) * (ky(2) - ky(1));

SpecE = cos(kphi) .* ExFun(kx, ky, kphi, kr) ...
    + sin(kphi) .* EyFun(kx, ky, kphi, kr);

%% Plotting
figure;
showImage(kx, ky, SpecE, DisplayFormat="dB", ColorScaleRangeDB=80);






