function mustBeFractional(val)
%Verify that value is between 0 and 1, exclusive.
% This is for nLayerViewer.
% 
% Author: Matt Dvorsky

mustBeInRange(val, 0, 1, "exclusive");

end