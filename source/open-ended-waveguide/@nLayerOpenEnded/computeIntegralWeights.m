function [] = computeIntegralWeights(self)
%Compute weights and nodes for mode integrals.
% This function is called whenever a parameter changes that would change
% the mode spectrums.
%
% Author: Matt Dvorsky

arguments
    self nLayerOpenEnded;
end

%% Fixed Point Integration Weights and Nodes
[self.fixed_kr, self.fixed_Ah, self.fixed_Ae] = ...
    self.computeAhat();

self.fixed_Ah = reshape(self.fixed_Ah, numel(self.fixed_kr), 1, []);
self.fixed_Ae = reshape(self.fixed_Ae, numel(self.fixed_kr), 1, []);

%% Set Flag
self.shouldRecomputeWeights = false;

end





