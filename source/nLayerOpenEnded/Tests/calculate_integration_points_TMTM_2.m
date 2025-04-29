function [Jr, Jw] = calculate_integration_points_TMTM_2(modeStruct, N)
%CALCULATE_INTEGRATION_POINTS_TMTM Calculate radial function TM to TM.
%
% Author: Matt Dvorsky

%% Integration Points and Line Segments
[nodes(:, 1), weights(:, 1)] = fejer2(N, -0.5, 0.5);

xp1 = modeStruct.BoundaryPoints(:, 1);
yp1 = modeStruct.BoundaryPoints(:, 2);

xp2 = circshift(xp1, -1);
yp2 = circshift(yp1, -1);

angInt = atan2(yp2 - yp1, xp2 - xp1);

%% Calculate Values
AzN_fun = @(x, y, ang) sin(ang).*modeStruct.HertzAz_dx(x, y) ...
    - cos(ang).*modeStruct.HertzAz_dy(x, y);

for ii = 1:numel(xp1)
    xint{ii} = (0.5*(xp2(ii) + xp1(ii)) + nodes.*(xp2(ii) - xp1(ii)));
    yint{ii} = (0.5*(yp2(ii) + yp1(ii)) + nodes.*(yp2(ii) - yp1(ii)));

    weightsAll = weights .* hypot(xp2(ii) - xp1(ii), yp2(ii) - yp1(ii));
    AzN{ii} = weightsAll .* AzN_fun(xint{ii}, yint{ii}, angInt(ii));
end

%% Calculate Off-Diagonal Cross Products
for ii = 1:numel(xint)
    for jj = 1:numel(yint)
        Jr{ii, jj} = reshape(hypot(xint{ii} - xint{jj}.', yint{ii} - yint{jj}.'), [], 1);
        Jw{ii, jj} = (2*pi) .* reshape(AzN{ii} .* AzN{jj}.', [], 1);
    end
end

%% Calculate Diagonal Cross Products
nodesY = 0.5*(0.5 + nodes.') + 0*nodes;
nodesX = (1 - 2*nodesY) .* nodes;

weightsAll = (weights .* weights.') .* (1 - 2*nodesY);

nodes1 = (nodesX + nodesY);
nodes2 = (nodesX - nodesY);

nodes1 = [nodes1(:); -nodes1(:)];
nodes2 = [nodes2(:); -nodes2(:)];
weightsAll = [weightsAll(:); weightsAll(:)];

% nodes1 = nodes + 0*nodes.';
% nodes2 = 0*nodes + nodes.';
% weightsAll = weights .* weights.';
% nodes1 = nodes1(:);
% nodes2 = nodes2(:);
% weightsAll = weightsAll(:);








% figure;
% plot(nodes1, nodes2, ".", LineWidth=1.5);
% grid on;
% axis image;
% 
% figure;
% plot3(nodes1, nodes2, weightsAll, ".", LineWidth=1.5);
% grid on;
% axis image;


% figure;
% plot(Jw{1, 1});
% tmp = Jw{1, 1};


for ii = 1:numel(xp1)
    xc = 0.5*(xp1(ii) + xp2(ii));
    yc = 0.5*(yp1(ii) + yp2(ii));
    xlen = xp2(ii) - xp1(ii);
    ylen = yp2(ii) - yp1(ii);

    x1 = xc + xlen.*nodes1;
    y1 = yc + ylen.*nodes1;
    x2 = xc + xlen.*nodes2;
    y2 = yc + ylen.*nodes2;

    Az1 = AzN_fun(x1, y1, angInt(ii));
    Az2 = AzN_fun(x2, y2, angInt(ii));

    Jr{ii, ii} = hypot(x1 - x2, y1 - y2);
    Jw{ii, ii} = (2*pi) .* weightsAll .* (xlen.^2 + ylen.^2) .* (Az1 .* Az2);
end


% figure;
% plot(Jw{1, 1});

%% Format Outputs
Jr = cell2mat(Jr(:));
Jw = cell2mat(Jw(:));

end

