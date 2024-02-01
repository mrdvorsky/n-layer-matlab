function [A] = computeA(O, f, er, ur, thk)
%COMPUTEA Compute the matrix A for each frequency.
% This function computes the matrix A as a function of each frequency
% specified by "f", which is used to compute the unnormalized mode
% S-parameter matrix.
%
% Example Usage:
%   [A] = O.computeA(f, er, ur, thk);
%   [K] = O.computeK(f);
%   idMat = eye(size(A, 1));
%   Smn = pagemldivide(idMat + A.*K, idMat - A.*K);
%
% Although the example above shows usage with a scalar value for "f", the
% input "f" can be a vector. In this case, the size of the 3rd dimension
% of each output matrix will be equal to numel(f).
%
% Inputs:
%
% Outputs:
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

%% Format Output
A = reshape(A, O.numModes, O.numModes, length(k0));

end

