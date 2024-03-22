function [] = computeIntegralWeights(O)
%COMPUTEINTEGRALWEIGHTS Compute weights and nodes for integral.
% This function is called whenever a parameter changes that would change
% the mode spectrums.
%
% Author: Matt Dvorsky

arguments
    O;
end

%% Fixed Point Integration Weights and Nodes
[O.fixed_kr, O.fixed_Ah, O.fixed_Ae] = ...
    O.computeAhat();

O.fixed_Ah = reshape(O.fixed_Ah, numel(O.fixed_kr), 1, []);
O.fixed_Ae = reshape(O.fixed_Ae, numel(O.fixed_kr), 1, []);

%% Set Flag
O.shouldRecomputeWeights = false;

end





