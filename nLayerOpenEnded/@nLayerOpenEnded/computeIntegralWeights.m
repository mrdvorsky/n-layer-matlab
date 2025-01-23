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

%% SVD Approximation
h_tmp = reshape(O.fixed_Ah(:, :).', O.numModes, O.numModes, []);
e_tmp = reshape(O.fixed_Ae(:, :).', O.numModes, O.numModes, []);

[Uh, Sh, Vh] = pagesvd(h_tmp);
[Ue, Se, Ve] = pagesvd(e_tmp);

Sh = Sh .* diag((1:O.numModes) <= 4);
Se = Se .* diag((1:O.numModes) <= 4);

h_new = pagemtimes(Uh, pagemtimes(Sh, pagectranspose(Vh)));
e_new = pagemtimes(Ue, pagemtimes(Se, pagectranspose(Ve)));
err_h = h_tmp - h_new;
err_e = e_tmp - e_new;

fixed_Ah_new = reshape(permute(h_new, [3, 1, 2]), size(O.fixed_Ah));
fixed_Ae_new = reshape(permute(e_new, [3, 1, 2]), size(O.fixed_Ae));

% O.fixed_Ah = fixed_Ah_new;
% O.fixed_Ae = fixed_Ae_new;

%% Set Flag
O.shouldRecomputeWeights = false;

end





