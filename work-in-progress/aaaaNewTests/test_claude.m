clc; clear; close all;

%% Inputs
NL = nLayerCircular(0, 1, waveguideR=0.5*4.81, modeSymmetryAxial="TM");
f = 300;


kc = NL.mode_kc0;
wgR = NL.waveguideR;
k_val = 1.5 - 0.001j;

sqrt_bc = @(kr) sqrt_with_branch(kr, k_val);

function kz = sqrt_with_branch(kr, k)
    kz = sqrt(k^2 - kr.^2);
    if imag(kz) > 0
        kz = -kz;
    end
end

%% The integrand near the pole
% f(kr) = g(kr) * J0(...) * H0(...) * kr^3 / (kr^2 - kc^2)^2
% At kr = kc, we have a second-order pole

% For second-order pole at kr = kc:
% Residue = lim_{kr->kc} d/dkr[(kr-kc)^2 * f(kr)]

% Let's compute this symbolically/numerically
syms kr_sym
assume(kr_sym, 'real');

% The part that's smooth at kr = kc
smooth_part = @(kr) sqrt_bc(kr) .* besselj(0, wgR*kr) .* ...
                    kr.^3 ./ (kr + kc).^2;

% For H0^(1) term
smooth_part_H1 = @(kr) sqrt_bc(kr) .* besselj(0, wgR*kr) .* ...
                       besselh(0, 1, wgR*kr) .* kr.^3 ./ (kr + kc).^2;

% For H0^(2) term  
smooth_part_H2 = @(kr) sqrt_bc(kr) .* besselj(0, wgR*kr) .* ...
                       besselh(0, 2, wgR*kr) .* kr.^3 ./ (kr + kc).^2;

% Residue = derivative of smooth part at kr = kc
h = 1e-7;
residue_H1 = (smooth_part_H1(kc + h) - smooth_part_H1(kc - h)) / (2*h);
residue_H2 = (smooth_part_H2(kc + h) - smooth_part_H2(kc - h)) / (2*h);

fprintf('Residue for H0^(1) term: %.8f + %.8fi\n', real(residue_H1), imag(residue_H1));
fprintf('Residue for H0^(2) term: %.8f + %.8fi\n', real(residue_H2), imag(residue_H2));
fprintf('Sum of residues:         %.8f + %.8fi\n', real(residue_H1 + residue_H2), imag(residue_H1 + residue_H2));

%% Now verify with contour integration
radius = 0.01;
circle_H1 = @(theta) smooth_part_H1(kc + radius*exp(1j*theta)) ./ ...
                     (kc + radius*exp(1j*theta) - kc).^2 .* ...
                     (1j*radius*exp(1j*theta));
                     
circle_H2 = @(theta) smooth_part_H2(kc + radius*exp(1j*theta)) ./ ...
                     (kc + radius*exp(1j*theta) - kc).^2 .* ...
                     (1j*radius*exp(1j*theta));

contour_res_H1 = integral(circle_H1, 0, 2*pi, 'ArrayValued', true) / (2*pi*1j);
contour_res_H2 = integral(circle_H2, 0, 2*pi, 'ArrayValued', true) / (2*pi*1j);

fprintf('\nContour residue H0^(1):  %.8f + %.8fi\n', real(contour_res_H1), imag(contour_res_H1));
fprintf('Contour residue H0^(2):  %.8f + %.8fi\n', real(contour_res_H2), imag(contour_res_H2));

%% Now the integrals
g = @(kr) 1 ./ sqrt_bc(kr);
JH1 = @(x) besselj(0, x).^2;
JH_upper = @(x) 0.5 * besselj(0, x) .* besselh(0, 1, x);
JH_lower = @(x) 0.5 * besselj(0, x) .* besselh(0, 2, x);

pole_term = @(kr) (kr./(kr.^2 - kc.^2)).^2 .* kr;

int_real = integral(@(kr) g(kr) .* JH1(wgR*kr) .* pole_term(kr), 0, 5000);
int_upper = integral(@(kr) g(kr) .* JH_upper(wgR*kr) .* pole_term(kr), 0, 60*(1+1j));
int_lower = integral(@(kr) g(kr) .* JH_lower(wgR*kr) .* pole_term(kr), 0, 60*(1-1j));

fprintf('\n--- Integrals ---\n');
fprintf('Real axis:  %.6f + %.6fi\n', real(int_real), imag(int_real));
fprintf('Upper:      %.6f + %.6fi\n', real(int_upper), imag(int_upper));
fprintf('Lower:      %.6f + %.6fi\n', real(int_lower), imag(int_lower));
fprintf('Sum:        %.6f + %.6fi\n', real(int_upper + int_lower), imag(int_upper + int_lower));
fprintf('Difference: %.6f + %.6fi\n', real(int_real - int_upper - int_lower), imag(int_real - int_upper - int_lower));

%% What should the residue contribution be?
% If upper contour encloses pole CCW: +2*pi*i*res
% If lower contour encloses pole CW: -2*pi*i*res (or +2*pi*i*res depending on orientation)

% Let's test:
test1 = int_upper + int_lower + 2*pi*1j*(residue_H1 + residue_H2);
test2 = int_upper + int_lower - 2*pi*1j*(residue_H1 + residue_H2);
test3 = int_upper + int_lower + 2*pi*1j*residue_H1 - 2*pi*1j*residue_H2;
test4 = int_upper + int_lower - 2*pi*1j*residue_H1 + 2*pi*1j*residue_H2;

fprintf('\n--- Testing residue signs ---\n');
fprintf('Real axis:           %.6f + %.6fi\n', real(int_real), imag(int_real));
fprintf('Test (+both):        %.6f + %.6fi\n', real(test1), imag(test1));
fprintf('Test (-both):        %.6f + %.6fi\n', real(test2), imag(test2));
fprintf('Test (+H1, -H2):     %.6f + %.6fi\n', real(test3), imag(test3));
fprintf('Test (-H1, +H2):     %.6f + %.6fi\n', real(test4), imag(test4));