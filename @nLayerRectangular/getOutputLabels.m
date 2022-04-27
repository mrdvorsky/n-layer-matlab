function [outputLabels] = getOutputLabels(O)
%GETOUTPUTLABELS Return a list of output channel labels.
% If gam is the output of "calculateGamma", outputLabels(ii) describes the
% set of measurements gam(:, ii).
%
% Example Usage:
%   NL = nLayerRectangular(maxM, maxN, band=wgBand);
%   outputLabels = NL.getOutputLabels();
%
% Outputs:
%   outputLabels - Vector of strings labeling each output channel.
%       Currently only outputs ["S11"].
%
% Author: Matt Dvorsky

arguments
    O;
end

%% Return Output Channel List
outputLabels = "S11";

end

