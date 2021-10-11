function [A1b1_EH] = integralFunction(O, tauP, k0, er, ur, thk)
%CONSTRUCTMATRIXEQUATION Summary of this function goes here
%   Detailed explanation goes here

tau = O.integralScaleFactor * (1 - tauP) ./ tauP;
[specE, specH] = O.multilayerSpectrumRect(tau, k0, er, ur, thk);

% if length(tau) == size(O.init_A1b1_E, 1)
%     A1b1_EH = specE .* O.init_A1b1_E + specH .* O.init_A1b1_H;
%     return;
% end

% Get indices and mixing factors for linear interolation
fracInd = tauP * (O.interpolationPointsTau - 1) + 1;
intInd = floor(fracInd);
m = fracInd - intInd;

% Perform linear interpolation
vLower = O.A1b1_EH(intInd, :, :, :);
vHigher = O.A1b1_EH(intInd + 1, :, :, :);
interp_A1b1_EH = vLower + m .* (vHigher - vLower);

A1b1_EH = specE.*interp_A1b1_EH(:, :, :, 1) ...
    + specH.*interp_A1b1_EH(:, :, :, 2);

end


