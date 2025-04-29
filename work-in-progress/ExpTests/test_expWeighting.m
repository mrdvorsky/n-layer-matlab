clc;
clear;
close all;

%% Inputs
f = 35;
er = {4 - 0.5j};
ur = {1};
thk = {1.5};

%% nLayer
NL = nLayerCircular(0, 10, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");
NL.receiveModeIndices = 1:NL.numModes;
gam = NL.calculate(f, er, ur, thk);

k0 = 2*pi * f ./ NL.speedOfLight;
beta = sqrt(k0.^2 - NL.mode_kc0.^2);

%% Plot
kr(:, 1) = linspace(0, 10, 1001);

spec1 = NL.modeStructs(1).EySpec(0, 0, kr, 0);

Ey1 = 0;
for ii = NL.receiveModeIndices
    specFun = NL.modeStructs(ii).EySpec;
    Ey1 = Ey1 + (1 - gam(ii)) .* specFun(0, 0, kr, 0) .* spec1;
end

Ey2 = 0;
[Gh, ~] = nLayer.computeGamma0(kr, 2*pi*f ./ NL.speedOfLight, er, ur, thk);
for ii = NL.receiveModeIndices
    specFun = NL.modeStructs(ii).EySpec;
    Ey2 = Ey2 + (k0 ./ beta(ii)) * (1 + gam(ii)) .* specFun(0, 0, kr, 0) .* spec1 .* Gh;
end

[Gh, ~] = nLayer.computeGamma0(kr, 2*pi*f ./ NL.speedOfLight, er, ur, thk);
for ii = NL.receiveModeIndices
    specFun = NL.modeStructs(ii).EySpec;
    EyFit(:, ii) = (1 + gam(ii)) .* specFun(0, 0, kr, 0) .* spec1 .* Gh;
end

x = EyFit \ Ey1;
Ey2 = EyFit * x;


figure;
plot(kr, real(Ey1), "", LineWidth=1.5);
hold on;
plot(kr, real(Ey2), ":", LineWidth=1.5);
grid on;

figure;
plot(kr, imag(Ey1), "", LineWidth=1.5);
hold on;
plot(kr, imag(Ey2), ":", LineWidth=1.5);
grid on;





