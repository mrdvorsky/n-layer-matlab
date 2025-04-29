clc;
clear;
close all;

%% Inputs
NL = nLayerRectangular(3, 2, waveguideBand="Ka");
modeStruct = NL.modeStructs(end);
kc = modeStruct.CutoffWavenumber;

ExSpec = modeStruct.ExSpec;
EySpec = modeStruct.EySpec;
W_TM = @(kr, kphi) (cos(kphi).*ExSpec(kr.*cos(kphi), kr.*sin(kphi), 0, 0) ...
    + sin(kphi).*EySpec(kr.*cos(kphi), kr.*sin(kphi), 0, 0)).^2;

%% Integration 1
[kphi(1, 1, :), kphi_w(1, 1, :)] = fejer2(128, 0, pi/2);
kphi_w = 4 * kphi_w;

W_kr1 = @(kr) sum(W_TM(kr, kphi) .* kphi_w, 3);

int1 = integral(W_kr1, 0, inf);

%% Integration 2
[Jr(1, 1, :), Jw(1, 1, :)] = calculate_integration_points_TMTM_2(modeStruct, 32);
W_kr2 = @(kr) sum(Jw .* besselj(0, Jr .* kr), 3) ...
    .* kr.^2 ./ (kr.^2 - kc.^2).^2;
% W_kr2 = @(kr) sum(Jw .* besselj(0, Jr .* kr), 3);

%% Actual Integrals
int1 = integral(W_kr1, 10 + 1j, 11, RelTol=1e-8)
int2 = integral(W_kr2, 10 + 1j, 11, RelTol=1e-8)

err_db = db((int1 - int2) ./ int1)

%% Plotting
% krPlot = linspace(0.1, 200, 200);

% spec1 = W_kr1(krPlot);
% spec2 = W_kr2(krPlot);

% figure;
% plot(krPlot, spec1, "", LineWidth=1.5);
% grid on;

% figure;
% plot(krPlot, spec2, "", LineWidth=1.5);
% grid on;
% 
% figure;
% plot(krPlot, db(spec1 - spec2), "", LineWidth=1.5);
% grid on;
% ylim([-300, 0]);














