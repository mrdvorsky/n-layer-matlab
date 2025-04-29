function [x0, coeff_k] = polyFourierTransform(x1, x2, vals)
%POLYFOURIERTRANSFORM Analytic Fourier transform of piecewise polynomials.
%
%   spec = sum(nufft(c, x0 ./ (2*pi), kx, 1) ...
%       .* (1./kx).^(1:size(c, 2)), 2);
%
% Author: Matt Dvorsky

arguments
    x1(:, 1);
    x2(:, 1) {mustHaveEqualSizes(x1, x2)};
    vals(:, :) {mustHaveEqualSizes(vals, x2, Dimensions=1)};
end

%% Calculate Change of Basis
lagPolys = lagrangePolys(size(vals, 2) - 1);
[fourierPolysA, fourierPolysB] = fourierPolys(size(vals, 2) - 1);

basisA = (fourierPolysA * lagPolys).';
basisB = (fourierPolysB * lagPolys).';

%% Compute Coefficients
scaleValues = 1 ./ (0.5*(x2 - x1)).^(0:size(vals, 2) - 1);
scaleValues(abs(x2 - x1) < 1e-10, 2:end) = 0;

coeff_k = [scaleValues .* (vals * basisA); ...
    scaleValues .* (vals * basisB)];
x0 = [x1; x2];

end



%% Helper Functions
function [lagPolys] = lagrangePolys(orderNum)
    n_plus1 = orderNum + 1;

    lagPolys = zeros(n_plus1, n_plus1);
    lagPolys(1, :) = 1;
    lagPolys(2, :) = linspace(1, -1, n_plus1);
    
    legPolysSpec = fft(lagPolys, [], 1);
    for ii = 1:n_plus1
        inds = [1:ii - 1, ii + 1:n_plus1];
        scale = 0.5*(n_plus1 - 1) ./ [ii - 1:-1:1, -(1:n_plus1 - ii)];
    
        lagPolys(:, ii) = prod(scale .* legPolysSpec(:, inds), 2);
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



