function [specQuery] = extrapolateSpectrum(specFun, Lt, krMax, N, krQuery)
%EXTRAPOLATESPECTRUM Extrapolate spectrum function.
% 
% Author: Matt Dvorsky

arguments
    specFun(1, 1);
    Lt(1, 1);
    krMax(1, 1);
    N(1, 1);
    krQuery(:, 1);
end

%% Sample Function
krSamp(:, 1) = linspace(-krMax, krMax, 2*N - 1);
rSamp = fftCoordinates(krSamp, ApplyFftShift=true);

% specSamp = specFun(krSamp) ./ sinc(0.5 * krSamp ./ krMax).^2;
specSamp = specFun(krSamp);
spatSamp = fftshift(ifft(ifftshift(specSamp)));

%% Query
specQuery = nufft(spatSamp, rSamp ./ (2*pi), krQuery) .* sinc(0.5 * krQuery ./ krMax).^2;



end

