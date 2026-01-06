clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(0, 1, waveguideR=0.5*4.81, modeSymmetryAxial="TM");
f = 300;

%% Modes
wgR = NL.waveguideR;        % Value is 2.41
kc = NL.mode_kc0;           % Value is 1

%% Integrands
JH1 = @(x) besselj(0, x).^2;
JH1c = @(x) besselj(0, x, 1).^2;
JH2 = @(x) 0.25*(besselh(0, 1, x) + besselh(0, 2, x)).^2;

intFun1 = @(kr) (JH1(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;
intFun2 = @(kr) (JH2(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;


JH2_1 = @(x) 0.25*(besselh(0, 1, x)).^2;
JH2_2 = @(x) 0.25*(besselh(0, 2, x)).^2;
JH2_3 = @(x) 0.5*(besselh(0, 1, x).*besselh(0, 2, x));

intFun2_1 = @(kr) (JH2_1(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;
intFun2_2 = @(kr) (JH2_2(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;
intFun2_3 = @(kr) (JH2_3(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;


%% Test
kMax = 5;
intFun5 = @(kr) (JH1c(wgR.*kr) .* exp(-2*kMax)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;

%% Full Integral
sqrt_bc = @(x) sqrt_continuous_vec(x);
g = @(kr) 1 ./ sqrt_bc((2.5 - 0.1j).^2 - kr.^2);

int1 = integral(@(kr) g(kr) .* intFun1(kr), 0, 50000)
int2 = integral(@(kr) g(kr) .* intFun2(kr), 0, 60*(1 + 0j))


int2_1 = integral(@(kr) g(kr) .* intFun2_1(kr), 0, 60*(1 + 1j))
int2_2 = integral(@(kr) g(kr) .* intFun2_2(kr), 0, 60*(1 - 1j), Waypoints=[1.5 + 0.1j])
% int2_2 = integral(@(kr) g(kr) .* intFun2_2(kr), 0, 60*(1 - 1j))
int2_3 = integral(@(kr) g(kr) .* intFun2_3(kr), 0, 60*(1 + 1j))
% int2_3 = integral(@(kr) g(kr) .* intFun2_3(kr), 0, 60*(1 + 0j), Waypoints=[1 - 1j])

int5_1 = integral(@(kr) g(kr) .* intFun5(kr), 0, kMax*(1 + 1j)) * exp(2*kMax)
int5_2 = integral(@(kr) g(kr) .* intFun5(kr), kMax*(1 + 1j), 20) * exp(2*kMax)


%% Residue
% int_res2 = integral(@(kr) g(kr) .* intFun2(kr), 0.5, 0.5, Waypoints=[1.1 + 1j, 1.1 - 1j])
% int_res3 = integral(@(kr) g(kr) .* intFun2c(kr), 0.5, 0.5, Waypoints=[1.1 + 1j, 1.1 - 1j])

% int_res2_1 = integral(@(kr) g(kr) .* intFun2_1(kr), 0.5, 0.5, Waypoints=[1.1 + 1j, 1.1 - 1j])
% int_res2_2 = integral(@(kr) g(kr) .* intFun2_2(kr), 0.5, 0.5, Waypoints=[1.1 + 1j, 1.1 - 1j])
% int_res2_3 = integral(@(kr) g(kr) .* intFun2_3(kr), 0.5, 0.5, Waypoints=[1.1 + 1j, 1.1 - 1j])


%% Helper
% function val = sqrt_continuous_vec(x)
%     val = sqrt(x);
%     % Flip sign where imaginary part is positive
%     mask = imag(val) > 0;
%     val(mask) = -val(mask);
% end

function val = sqrt_continuous_vec(x)
    val = sqrt(x);
    val = real(val) - 1j*abs(imag(val));
end