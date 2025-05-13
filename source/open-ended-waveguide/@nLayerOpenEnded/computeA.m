function [A] = computeA(self, f, er, ur, thk)
%Compute the matrix A for each frequency.
% This function computes the matrix A as a function of each frequency
% specified by "f", which is used to compute the mode S-parameter matrix.
%
% Example Usage:
%   [A] = self.computeA(f, er, ur, thk);
%   [K] = self.computeK(f);
%   idMat = eye(size(A, 1));
%   Smn = pagemldivide(idMat + A.*K, idMat - A.*K);
%
% Author: Matt Dvorsky

arguments
    self nLayerOpenEnded;

    f;
    er;
    ur;
    thk;
end

%% Calculate A
k0(1, 1, 1, :) = 2*pi .* f(:) ./ self.speedOfLight;
[Gamma0h, Gamma0e] = nLayer.computeGamma0(self.fixed_kr, k0, er, ur, thk);

% Gamma0h = Gamma0h .* (0.5 - 0.5*tanh(O.fixed_kr - 10)).^2;
% Gamma0e = Gamma0e .* (0.5 - 0.5*tanh(O.fixed_kr - 10)).^2;
A = pagemtimes(Gamma0h, "transpose", self.fixed_Ah, "none") ...
    + pagemtimes(Gamma0e, "transpose", self.fixed_Ae, "none");

% figure;
% plots(real(Gamma0h), "", LineWidth=1.5);
% hold on;
% plots(imag(Gamma0h), "", LineWidth=1.5);

%% Format Output
A = reshape(A, self.numModes, self.numModes, numel(k0));

end

