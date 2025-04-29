clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(0, 1, waveguideBand="Ka_TE01", modeSymmetryAxial="TM");
modeStruct = NL.modeStructs(end);
kc = modeStruct.CutoffWavenumber;

kx(:, 1) = linspace(-0.5, 5, 600);
ky(1, :) = linspace(-1.5, 1.5, 200);

%% Calculate 
modeFun = @(x) modeStruct.ExSpec(0, 0, x, 0);

GamH = modeFun(kx + 1j*ky);

%% Plot
figure;
ImgH = showImage(kx, ky, GamH, DisplayFormat="Magnitude");
grid on;
clim(10*[0, 1]);

figure;
ImgH = showImage(kx, ky, GamH, DisplayFormat="Phase");
grid on;
% clim(10*[0, 1]);






