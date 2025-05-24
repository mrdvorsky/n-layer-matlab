function [er, ur, thk] = arrayToStructure(structureArray)
%Convert structure from 2D array form to "structure form".
% This is for nLayerViewer.
% 
% Author: Matt Dvorsky

arguments
    structureArray(:, :);
end

%% Convert
er = num2cell(structureArray(:, 1) .* (1 - 1j*structureArray(:, 2)));
ur = num2cell(structureArray(:, 3) .* (1 - 1j*structureArray(:, 4)));
thk = num2cell(structureArray(:, 5));

end
