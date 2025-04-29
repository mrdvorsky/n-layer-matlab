function [x, y, weights] = gaussCircle(rOuter, numR, numPhi, rInner)
%GAUSSCIRCLE Summary of this function goes here
%   Detailed explanation goes here

arguments
    rOuter;
    numR;
    numPhi;
    rInner = 0;
end


[r(:, 1), weights_r(:, 1)] = fejer2(numR, rInner, rOuter);
[phi(1, :), weights_phi(1, :)] = fejer2(numPhi, 0, 0.5*pi);
phi = [phi, phi + 0.5*pi, phi + 1.0*pi, phi + 1.5*pi];
weights_phi = repmat(weights_phi, 1, 4);

x = r .* cos(phi);
y = r .* sin(phi);
weights = r .* weights_r .* weights_phi;

x = x(:);
y = y(:);
weights = weights(:);

end

