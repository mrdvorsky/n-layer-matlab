clc;
clear;
close all;

%% Inputs
er = [2.1 - 0.001j];
ur = [1];
thk = [10];

f = linspace(8.2, 12.4, 801);

%% Integral
% intVal1 = integral(@(k) S_int(k, 4, 1).^2, -inf, inf)
% intVal2 = integral(@(k) C_int(k, 4, 2).^2, -inf, inf)
% return;

%% nLayerRectangular
tic;
NL1 = nLayerRectangularOld(3, 2, waveguideBand="x", convergenceAbsTol=1e-4, verbosity=1);
toc;


%% nLayerGeneral
wgA = NL1.waveguideA/1;
wgB = NL1.waveguideB/1;

NL2 = nLayerOpenEnded(verbosity=1, convergenceAbsTol=1e-4);
% NL2.waveguideEr = 4 - 0.4j;
% NL2.waveguideUr = 4 - 0.4j;

modesTE = NL1.modesTE;
for ii = 1:size(modesTE, 1)
    m = modesTE(ii, 1);
    n = modesTE(ii, 2);

    betaC = hypot(m*pi/wgA, n*pi/wgB);
    NL2.modeBetaCutoffTE(ii) = betaC;

    scaleFactor = pi / (betaC * sqrt(wgA*wgB)) ...
        ./ sqrt(1 + 0*(n == 0));

    modeSpectrumEx_TE{ii} = @(kx, ky, ~, ~) -(n/wgB) * scaleFactor ...
        .* C_int(kx, wgA, m) .* S_int(ky, wgB, n);

    modeSpectrumEy_TE{ii} = @(kx, ky, ~, ~) (m/wgA) * scaleFactor ...
        .* S_int(kx, wgA, m) .* C_int(ky, wgB, n);
end

modesTM = NL1.modesTM;
for ii = 1:size(modesTM, 1)
    m = modesTM(ii, 1);
    n = modesTM(ii, 2);

    betaC = hypot(m*pi/wgA, n*pi/wgB);
    NL2.modeBetaCutoffTM(ii) = betaC;

    scaleFactor = pi / (betaC * sqrt(wgA*wgB)) ...
        ./ sqrt(1 + 0*(n == 0));

    modeSpectrumEx_TM{ii} = @(kx, ky, ~, ~) (m/wgA) * scaleFactor ...
        .* C_int(kx, wgA, m) .* S_int(ky, wgB, n);

    modeSpectrumEy_TM{ii} = @(kx, ky, ~, ~) (n/wgB) * scaleFactor ...
        .* S_int(kx, wgA, m) .* C_int(ky, wgB, n);
end



tic;
NL2.recomputeInterpolants(...
    SpecEx_TE=modeSpectrumEx_TE, ...
    SpecEy_TE=modeSpectrumEy_TE, ...
    SpecEx_TM=modeSpectrumEx_TM, ...
    SpecEy_TM=modeSpectrumEy_TM, ...
    ModeSymmetryX="Even", ModeSymmetryY="Even", ...
    IntegralScaleFactor=pi.^2 ./ wgA);
toc;

%% Calculate
gam1 = NL1.calculate(f, er, ur, thk);
gam2 = NL2.calculate(f, er, ur, thk);

relErr = abs(max(gam1 - gam2))

%% Plot
figure;
plot(gam1, "-", LineWidth=1.5);
hold on;
plot(gam2, "-", LineWidth=1.5);
zplane([]);
grid on;
legend(["Original", "New"]);

% figure;
% plot(f, real(gam2), "-", LineWidth=1.5);
% hold on;
% plot(f, imag(gam2), "-", LineWidth=1.5);
% grid on;
% legend(["Real", "Imag"]);




%% Helper Functions
function v = S_int(k, a, m)

if m == 0
    v = 0;
    return;
end
v = sqrt( pi ) * m .* sinc((0.5/pi) .* (a.*k - pi*m)) ./ (k  + m.*pi./a);

end

function v = C_int(k, a, m)

v = sqrt(1/pi)*a * k .* sinc((0.5/pi) .* (a.*k - pi*m)) ./ (k  + m.*pi./a);
if m == 0
    v = v .* sqrt(0.5);
end

end

