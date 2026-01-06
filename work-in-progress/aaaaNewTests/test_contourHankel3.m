clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(0, 1, waveguideR=0.5*4.81, modeSymmetryAxial="TM");
f = 300;

%% Modes
specE1 = NL.waveguideModes.WeSpec;

wgR = NL.waveguideR;        % Value is 2.41
kc = NL.mode_kc0;           % Value is 1

%% Integrands
% JH1 = @(x) besselj(0, x, 1).*besselh(0, 1, x, 1).*exp(1j*real(x));
% JH2 = @(x) besselj(0, x, 1).*besselh(0, 2, x, 1).*exp(-1j*real(x));
JH1 = @(x) besselj(0, x).^2;
JH2 = @(x) 0.5*besselj(0, x).*besselh(0, 1, x);
JH2c = @(x) 0.5*besselj(0, x).*besselh(0, 2, x);

intFun1 = @(kr) (JH1(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;
intFun2 = @(kr) (JH2(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;
intFun2c = @(kr) (JH2c(wgR.*kr)).*(kr./(kr.^2 - kc.^2)).^2 .* kr;

%% Full Integral
% sqrt_bc = @(x) -conj(sqrt(x)) .* sign(imag(sqrt(x)) + eps);
sqrt_bc = @(x) sqrt_continuous_vec(x);
g = @(kr) 1 ./ sqrt_bc((1.5 - 0.001j).^2 - kr.^2);

int1 = integral(@(kr) g(kr) .* intFun1(kr), 0, 50000)
int2 = integral(@(kr) g(kr) .* intFun2(kr), 0, 60*(1 + 1j))
int3 = integral(@(kr) g(kr) .* intFun2c(kr), 0, 500 + 0.1j)


int_res2 = integral(@(kr) g(kr) .* intFun2(kr), 0.5, 0.5, Waypoints=[1.1 + 1j, 1.1 - 1j])
int_res3 = integral(@(kr) g(kr) .* intFun2c(kr), 0.5, 0.5, Waypoints=[1.1 + 1j, 1.1 - 1j])



%%
x(:, 1) = 5*linspace(-2, 2, 300);
y(1, :) = 5*linspace(-2, 2, 300);

figure;
showImage(x, y, g(x + 1j*y));

%%



% function val = sqrt_continuous_vec(x)
%     val = exp(1j*deg2rad(45)).*sqrt(exp(-1j*deg2rad(90))*x);
%     % Flip sign where imaginary part is positive
%     % mask = imag(val) > 0;
%     % val(mask) = -val(mask);
% end

function val = sqrt_continuous_vec(x)
    val = sqrt(x);
    % Flip sign where imaginary part is positive
    mask = imag(val) < 0;
    val(mask) = -val(mask);
end