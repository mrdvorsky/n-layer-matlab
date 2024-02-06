function [er, ur, thk] = changeStructureConductivity(O, f, er, ur, thk, sigma)
%CHANGESTRUCTURECONDUCTIVITY Modify multilayer structure with new backing conductivity.
% Wrapper around "nLayer.changeStructureConductivity", but the speed of
% light is passed in by default. See documentation of this function for
% more details.
%
% Author: Matt Dvorsky

arguments
    O;
    f(:, 1);
    er(:, :);
    ur(:, :);
    thk(1, :);
    sigma(:, 1);
end

%% Call "nLayer.changeStructureConductivity"
[er, ur, thk] = nLayer.changeStructureConductivity(f, er, ur, thk, sigma, ...
    SpeedOfLight=O.speedOfLight);

end
