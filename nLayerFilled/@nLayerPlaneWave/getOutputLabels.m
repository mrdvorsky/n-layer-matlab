function [outputLabels] = getOutputLabels(O)
%GETOUTPUTLABELS Return a list of output channel labels.
% If gam is the output of "calculate", outputLabels(ii) describes the
% set of measurements gam(:, ii).
%
% Example Usage:
%   NL = nLayerFilledRectangular(...);
%   outputLabels = NL.getOutputLabels();
%
% Outputs:
%   outputLabels - Vector of strings labeling each output channel.
%       Currently only outputs ["S_{11}", "S_{21}", "S_{12}", "S_{22}"].
%
% Author: Trent Moritz

arguments
    O;
end

%% Return Output Channel List
outputLabels = ["S_{11}", "S_{21}", "S_{12}", "S_{22}"];

if ~isempty(O.outputIndices)
    outputLabels = outputLabels(1, O.outputIndices);
end

end

