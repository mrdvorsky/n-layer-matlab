function [] = setInitialValues(O, options)
%SETINITIALVALUES Sets the initial values for each variable.
% This function is simply a setter for the initialValues_{...} parameters
% that checks array sizes. Any unspecified parameter indices will be left
% unchanged.
%
% Example Usage (for 2 layers):
%   NLsolver.setInitialValues(Er=[1 - 0.01j, 2 - 2j], Thk=[10, 20]);
%   NLsolver.setInitialValues(Ur=[2 - 1j, 1], Thk=[10, 1]);
%
% Author: Matt Dvorsky

arguments
    O;
    options.Er (1, :) {mustBeValidInitialValues(options.Er, O)}  = O.initialValue_er;
    options.Ur (1, :) {mustBeValidInitialValues(options.Ur, O)}  = O.initialValue_ur;
    options.Thk(1, :) {mustBeValidInitialValues(options.Thk, O)} = O.initialValue_thk;
end

%% Assign Initial Values
if iscell(options.Er)
    options.Er = cell2mat(options.Er);
end
if iscell(options.Ur)
    options.Ur = cell2mat(options.Ur);
end
if iscell(options.Thk)
    options.Thk = cell2mat(options.Thk);
end

O.initialValue_er  = options.Er;
O.initialValue_ur  = options.Ur;
O.initialValue_thk = options.Thk;

end

function mustBeValidInitialValues(vals, O)
    if numel(vals) ~= O.layerCount
        throwAsCaller(MException("nLayerInverse:mustBeValidInitialValues", ...
            "The parameters 'initialValue_{er, ur, thk}' must be " + ...
            "vectors with layerCount (%d) elements.", O.layerCount));
    end
end

