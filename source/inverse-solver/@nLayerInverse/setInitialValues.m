function [] = setInitialValues(self, options)
%Sets the initial values for each variable.
% This function is simply a setter for the initialValues_{...} parameters
% that checks array sizes. Any unspecified parameter indices will be left
% unchanged.
%
% Example Usage (for 2 layers):
%   NLsolver.setInitialValues(Er=[1 - 0.01j, 2 - 2j], Thk=[10, 20]);
%   NLsolver.setInitialValues(Ur=[2 - 1j, 1], Thk=[10, 1]);
%
%
% Author: Matt Dvorsky

arguments
    self nLayerInverse;

    options.Er (1, :) = self.initialValue_er;
    options.Ur (1, :) = self.initialValue_ur;
    options.Thk(1, :) = self.initialValue_thk;
end

%% Assign Initial Values
self.initialValue_er  = options.Er;
self.initialValue_ur  = options.Ur;
self.initialValue_thk = options.Thk;

end


