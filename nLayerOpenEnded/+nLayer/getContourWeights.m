function [krcNodes, krc, momentH_weights, momentE_weights] = getContourWeights(Nm, Nrho, L, Lc, Lch, Lcw)
%GETCONTOURWEIGHTS Helper functions for contour integration.
%
% Author: Matt Dvorsky

arguments
    Nm(1, 1);
    Nrho(1, 1);
    L(1, 1);
    Lc(1, 1);
    Lch(1, 1);
    Lcw(1, 1);
end

%% Helpers
krToKrc = @(x) x + 1j * (Lch.*Lcw) * x ./ sqrt((Lch + x.^2).*(3*Lcw.^2 + x.^2));
krToKrcPrime = @(x) 1 + 1j * (Lch.*Lcw) * (3*Lch*Lcw.^2 - x.^4) ./ ((Lch + x.^2).*(3*Lcw.^2 + x.^2)).^1.5;

%% Compute Weights and Nodes
% [kr(1, 1, 1, :), kr_weights(1, 1, 1, :)] = fejer2_halfOpen(Nrho, L);
% krc = krToKrc(kr);
% 
% moment_weights = krToKrcPrime(kr) .* kr_weights ...
%     .* sin(2*(1:Nm).' .* acot(sqrt(kr./Lc))) ...
%     ./ sin(2 * acot(sqrt(kr./Lc)));
% 
% momentH_weights = moment_weights ./ sqrt(1 + krc.^2);
% momentE_weights = moment_weights .* sqrt(1 + krc.^2);

[kr(1, 1, 1, :), kr_weights(1, 1, 1, :)] = fejer2_halfOpen(Nrho, L);
krc = kr;

krInv = conj(krToKrc(kr));
krInv = krInv - (krToKrc(krInv) - kr) ./ krToKrcPrime(krInv);
krInv = krInv - (krToKrc(krInv) - kr) ./ krToKrcPrime(krInv);
krInv = krInv - (krToKrc(krInv) - kr) ./ krToKrcPrime(krInv);
krInv = krInv - (krToKrc(krInv) - kr) ./ krToKrcPrime(krInv);
krInv = krInv - (krToKrc(krInv) - kr) ./ krToKrcPrime(krInv);
krInv = krInv - (krToKrc(krInv) - kr) ./ krToKrcPrime(krInv);
moment_weights = kr_weights ...
    .* sin(2*(1:Nm).' .* acot(sqrt(krInv./Lc))) ...
    ./ sin(2 *           acot(sqrt(krInv./Lc)));

momentH_weights = moment_weights ./ (1 + krc.^1);
momentE_weights = moment_weights .* (1 + krc.^1);

%% Compute Nodes
krNodes = fejer2_halfOpen(Nm, Lc);
krcNodes = krToKrc(krNodes);

end

