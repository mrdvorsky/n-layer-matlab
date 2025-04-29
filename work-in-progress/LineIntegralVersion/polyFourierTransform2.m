function [x0, s0, c0] = polyFourierTransform2(x1, x2, vals)
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
lagPolys = flip(lagrangePolys(size(vals, 2) - 1), 1);
legPolys = legendrePolynomials(size(vals, 2));

basisPoly = (legPolys.' \ lagPolys).';

%% Compute Coefficients
x0 = 0.5 * (x1 + x2);
s0 = 0.5 * (x2 - x1);
c0 = 2 * abs(s0) * (-1j).^(0:size(vals, 2) - 1) .* (vals * basisPoly);

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


