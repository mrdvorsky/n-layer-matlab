function [krContour, momentKernelFun] = nlGetContourIntegrals(N, L, Lc)
%GETCONTOURINTEGRALS Helper functions for contour integration.
%
% Author: Matt Dvorsky

arguments
    N(1, 1);
    L(1, 1);
    Lc(1, 1);
end

%% Helper Functions
% Functions to go to contour path (krc) and back (kr).
krToKrc = @(x) x + 1j .* Lc .* x ./ (Lc + (x));
krToKrcPrime = @(x) 1 + 1j .* Lc.^2 ./ ((x) + Lc).^2;

krcToKr = @(x) 0.5 * (sqrt(x.^2 + Lc*(2-2j)*x + 2j*Lc.^2) + x - Lc - 1j*Lc);
krcToKrPrime = @(x) 0.5 * (1 + (x + Lc - 1j*Lc) ...
    ./ sqrt(x.^2 + Lc*(2-2j)*x + 2j*Lc.^2));

%% Get Integral Nodes
[x(:, 1)] = nlFejer2_halfOpen(N, L);
krContour = krToKrc(x);

%% Define Kernel Functions for Moment Integrals
momentKernelFun = @(kr, n) krToKrcPrime(krcToKr(kr)) .* krcToKrPrime(kr) ...
    .* sin(2 .* n .* acot(sqrt(krcToKr(kr) ./ L))) ...
    ./ sin(2 .*      acot(sqrt(krcToKr(kr) ./ L)));

end

