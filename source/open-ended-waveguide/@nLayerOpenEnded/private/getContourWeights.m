function [krc, moment_weights] = getContourWeights(Nm, Nrho, a, b)
%GETCONTOURWEIGHTS Helper functions for contour integration.
%
% Author: Matt Dvorsky

arguments
    Nm(:, 1);
    Nrho(:, 1);
    a(:, 1);
    b(:, 1);
end

%% Compute Weights and Nodes for Each Interval
for ii = 1:numel(a)
    [x, w] = fejer2(Nrho{ii}, 0, pi);
    x = reshape(x, 1, 1, 1, []);
    w = reshape(w, 1, 1, 1, []);

    krc{ii} = 0.5*(b{ii} - a{ii}) .* cos(x) + 0.5*(a{ii} + b{ii});
    moment_weights{ii} = w .* sin(x .* (1:Nm{ii}).');
end

end

