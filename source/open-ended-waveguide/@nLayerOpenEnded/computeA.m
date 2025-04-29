function [A] = computeA(O, f, er, ur, thk)
%COMPUTEA Compute the matrix A for each frequency.
% This function computes the matrix A as a function of each frequency
% specified by "f", which is used to compute the mode S-parameter matrix.
%
% Example Usage:
%   [A] = O.computeA(f, er, ur, thk);
%   [K] = O.computeK(f);
%   idMat = eye(size(A, 1));
%   Smn = pagemldivide(idMat + A.*K, idMat - A.*K);
%
% Author: Matt Dvorsky

arguments
    O;
    f;
    er;
    ur;
    thk;
end

%% Calculate A
k0(1, 1, 1, :) = 2*pi .* f(:) ./ O.speedOfLight;
[Gamma0h, Gamma0e] = nLayer.computeGamma0(O.fixed_kr, k0, er, ur, thk);
A = pagemtimes(Gamma0h, "transpose", O.fixed_Ah, "none") ...
    + pagemtimes(Gamma0e, "transpose", O.fixed_Ae, "none");

% figure;
% plots(real(Gamma0h), "", LineWidth=1.5);
% hold on;
% plots(imag(Gamma0h), "", LineWidth=1.5);

%% Format Output
A = reshape(A, O.numModes, O.numModes, numel(k0));

end

