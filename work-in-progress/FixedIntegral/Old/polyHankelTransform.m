function [r0, coeff_k] = polyHankelTransform(r1, r2, vals, options)
%POLYFOURIERTRANSFORM Analytic Hankel transform of piecewise polynomials.
%
%   spec = sum(nufft(c, x0 ./ (2*pi), kx, 1) ...
%       .* (1./kx).^(1:size(c, 2)), 2);
%
% Author: Matt Dvorsky

arguments
    r1(:, 1) {mustBeNonnegative};
    r2(:, 1) {mustBeNonnegative, mustHaveEqualSizes(r1, r2)};
    vals(:, :) {mustHaveEqualSizes(vals, r2, Dimensions=1)};
    options.BesselOrder(1, 1) {mustBeReal} = 0;
end

%% Calculate Change of Basis
lagPolys = lagrangePolys(r1 ./ r2, size(vals, 2) - 1);
[fourierPolysA, fourierPolysB] = fourierPolys(size(vals, 2) - 1);

basisA = (fourierPolysA * lagPolys).';
basisB = (fourierPolysB * lagPolys).';

%% Compute Coefficients
scaleValues = 1 ./ (0.5*(r2 - r1)).^(0:size(vals, 2) - 1);

coeff_k = [scaleValues .* (vals * basisA); ...
    scaleValues .* (vals * basisB)];
r0 = [r1; r2];

end



%% Helper Functions
function [lagPolys] = lagrangePolys(ratio, orderNum)
    n_plus1 = orderNum + 1;
    xn = linspace(ratio, 1, n_plus1).^2;

    lagPolys = zeros(n_plus1, n_plus1);
    lagPolys(1, :) = 1;
    lagPolys(2, :) = -xn;
    
    legPolysSpec = fft(lagPolys, [], 1);
    for ii = 1:n_plus1
        inds = [1:ii - 1, ii + 1:n_plus1];
        scale = 1 ./ prod(xn(ii) - xn(inds));
    
        lagPolys(:, ii) = scale .* prod(legPolysSpec(:, inds), 2);
    end

    lagPolys = real(ifft(lagPolys, [], 1));
end

function [fourierPolysA, fourierPolysB] = fourierPolys(orderNum)
    fourierPolysA = zeros(orderNum + 1, orderNum + 1);
    for k = 0:orderNum
        fourierPolysA(1:k + 1, k + 1) = ...
            (-1j).^((k:-1:0) + k + 1) .* cumprod([1, k:-1:1]);
    end

    fourierPolysB = conj(fourierPolysA) .* (-1).^(0:orderNum);

    fourierPolysA = flip(fourierPolysA, 2);
    fourierPolysB = flip(fourierPolysB, 2);
end



