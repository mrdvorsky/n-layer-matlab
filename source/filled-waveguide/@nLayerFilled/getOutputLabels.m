function [outputLabels] = getOutputLabels(self)
%Return list of output labels for calculated reflection and transmission
% coefficients. 
% Since multimodal solutions are not supported, output labels will
% correspond to the Smn parameters calculated.
%
% Author: Trent Moritz
arguments
    self nLayerFilled
end

outputLabels = ["S_{11}", "S_{21}", "S_{12}", "S_{22}"];

end