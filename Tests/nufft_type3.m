function [y] = nufft_type3(x, t, f)
%NUFFT_TYPE3 Summary of this function goes here
%   Detailed explanation goes here



nt = numel(t);
dt = max(t) ./ (nt - 1);
ts = dt * (0:1:(nt - 1));

% fs = (1./dt) * (0:1:(nt - 1)) ./ (nt);
fs = fftCoordinates(1:nt, ApplyFftShift=true, PositiveOutput=false) * (0.5./pi);
% fs = fftCoordinates(ts, ApplyFftShift=false, PositiveOutput=true) * (0.5./pi);
t_scaled = t./dt;

tic;
% ys = ifftshift(nufft(x, t_scaled, fs));
ys = ifftshift(nufft(x .* exp(1j*pi*t_scaled), t_scaled, []));
% ys = nufft(x, t, fs);
toc;
tic;
xs = ifft(ys);
toc;
tic;
y = nufft(xs, ts, f) ./ numel(ts);
toc;


end

