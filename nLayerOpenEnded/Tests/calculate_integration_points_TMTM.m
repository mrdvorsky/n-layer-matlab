function [Jr, Jw] = calculate_integration_points_TMTM(modeStruct, N)
%CALCULATE_INTEGRATION_POINTS_TMTM Calculate radial function TM to TM.
%
% Author: Matt Dvorsky

%% Integration Points
[nodes(:, 1), weights(:, 1)] = fejer2(N, -0.5, 0.5);

xp(1, :) = modeStruct.BoundaryPoints(:, 1);
yp(1, :) = modeStruct.BoundaryPoints(:, 2);

xps = circshift(xp, -1);
yps = circshift(yp, -1);

xint = (0.5*(xps + xp) + nodes.*(xps - xp));
yint = (0.5*(yps + yp) + nodes.*(yps - yp));
angInt = atan2(yps - yp, xps - xp) + 0*xint;

weightsAll = weights .* hypot(xps - xp, yps - yp);

%% Calculate Values
xint = xint(:);
yint = yint(:);
angInt = angInt(:);
weightsAll = weightsAll(:);
AzN = sin(angInt).*modeStruct.HertzAz_dx(xint, yint) ...
    - cos(angInt).*modeStruct.HertzAz_dy(xint, yint);

Jr = reshape(hypot(xint - xint.', yint - yint.'), [], 1);
Jw = (2*pi) .* reshape((weightsAll.*AzN) .* (weightsAll.*AzN).', [], 1);

end

