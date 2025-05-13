function [outputLabels] = getOutputLabels(self)
%Return a list of output channel labels.
% If gam is the output of "calculateGamma", outputLabels(ii) describes the
% set of measurements gam(:, ii).
%
% Example Usage:
%   NL = nLayerRectangular(maxM, maxN, Band=wgBand);
%   outputLabels = NL.getOutputLabels();
%
%
% Outputs:
%   outputLabels - Vector of strings labeling each output channel.
%
% Author: Matt Dvorsky

arguments
    self nLayerOpenEnded;
end

%% Return Output Channel List
[modeLabelsRx, modeLabelsTx] = ndgrid(...
    self.modeLabels(self.receiveModeIndices), ...
    self.modeLabels(self.excitationModeIndices));

outputLabels = reshape(...
    compose("S_{%s,%s}", modeLabelsRx(:), modeLabelsTx(:)), ...
    [1, size(modeLabelsTx)]);

end

