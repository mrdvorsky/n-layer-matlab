function [krContour, momentKernelFun] = getContourIntegrals(N, L, Lc)
%GETCONTOURINTEGRALS Helper functions for contour integration.
%
% Author: Matt Dvorsky

arguments
    N(1, 1);
    L(1, 1);
    Lc(1, 1);
end

%% Helper Functions
% signReal = @(x) sign(real(x)) + (x == 0);

% Functions to go to contour path (krc) and back (kr).
krToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + (x));
krToKrcPrime = @(x) 1 + 1j .* Lc.^2 ./ ((x) + Lc).^2;

krcToKr = @(x) 0.5 * (sqrt(x.^2 + Lc*(2-2j)*x + 2j*Lc.^2) + x - Lc - 1j*Lc);
krcToKrPrime = @(x) 0.5 * (1 + (x + Lc - 1j*Lc) ...
    ./ sqrt(x.^2 + Lc*(2-2j)*x + 2j*Lc.^2));

% % Functions to go to contour path (krc) and back (kr).
% krToKrc = @(z) z .* exp(0.25j * pi ./ (1 + abs(z)));
% krToKrcPrime = @(z) exp(0.25j * pi ./ (1 + abs(z))) ...
%     .* (abs(z) + abs(z).^3 + (2 - 0.25j*pi)*abs(z).^2) ...
%     ./ abs(z) ./ (1 + abs(z)).^2;
% 
% krcToKr = @(z) z .* exp(-0.25j * pi ./ (1 + abs(z)));
% krcToKrPrime = @(z) exp(-0.25j * pi ./ (1 + abs(z))) ...
%     .* (abs(z) + abs(z).^3 + (2 + 0.25j*pi)*abs(z).^2) ...
%     ./ abs(z) ./ (1 + abs(z)).^2;

% % Functions to go to contour path (krc) and back (kr).
% krToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + (x));
% krToKrcPrime = @(x) 1 + 1j .* Lc.^2 ./ ((x) + Lc).^2;
% 
% krcToKr = @(x) 0.5 * (sqrt((x).^2 + Lc*(2-2j)*(x) + 2j*Lc.^2) + x - Lc - 1j*Lc);
% krcToKrPrime = @(x) 0.5 * (1 + ((x) + Lc - 1j*Lc) ...
%     ./ sqrt((x).^2 + Lc*(2-2j)*(x) + 2j*Lc.^2));

% % Alternate functions to go to contour path (krc) and back (kr).
% krToKrc = @(x) 0.5 * (sqrt(x.^2 + Lc*(2+2j)*x - 2j*Lc.^2) + x - Lc + 1j*Lc);
% krToKrcPrime = @(x) 0.5 * (1 + (x + Lc + 1j*Lc) ...
%     ./ sqrt(x.^2 + Lc*(2+2j)*x - 2j*Lc.^2));
% 
% krcToKr = @(x) x - 1j .* Lc .* x ./ (Lc + (x));
% krcToKrPrime = @(x) 1 - 1j .* Lc.^2 ./ ((x) + Lc).^2;

%% Get Integral Nodes
[x(:, 1)] = fejer2_halfOpen(N, L);
krContour = krToKrc(x);

%% Define Kernel Functions for Moment Integrals
momentKernelFun = @(kr, n) krToKrcPrime(krcToKr(kr)) .* krcToKrPrime(kr) ...
    .* sin(2 .* n .* acot(sqrt(krcToKr(kr) ./ L))) ...
    ./ sin(2 .*      acot(sqrt(krcToKr(kr) ./ L)));

end

