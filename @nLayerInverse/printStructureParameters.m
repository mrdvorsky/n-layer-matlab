function [varargout] = printStructureParameters(O, er, ur, thk, formatOptions, options)
%PRINTSTRUCTUREPARAMETERS Display multilayer structure solver parameters.
%   Detailed explanation goes here

arguments
    O;
    er(1, :) = [];
    ur(1, :) = [];
    thk(1, :) = [];
    
    formatOptions.ShowLimits(1, 1) {mustBeNumericOrLogical} = false;
    formatOptions.ShowInitialValues(1, 1) {mustBeNumericOrLogical} = true;
    formatOptions.SolveParameterFormatString {mustBeTextScalar} = "{%s}";
    
    options.Title {mustBeTextScalar} = "Solver Structure Parameters";
    options.BackingConductivity {mustBePositive} = inf;
    options.ThkUnitLabel {mustBeTextScalar} = "mm";
    options.ConductivityUnitLabel {mustBeTextScalar} = "S/m";
    
    options.Width(1, 1) {mustBeInteger, mustBePositive} = 72;
    options.ThkWidth(1, 1) {mustBeInRange(options.ThkWidth, 0, 1)} = 0.3;
    options.ErWidth(1, 1) {mustBeInRange(options.ErWidth, 0, 1)} = 0.4;
    options.UrWidth(1, 1) {mustBeInRange(options.UrWidth, 0, 1)} = 0.3;
    options.IndentWidth(1, 1) {mustBeInteger, mustBeNonnegative} = 4;
    options.FormatString(:, :) {mustBeText} = "%.4g";
    options.ConductivityFormatString {mustBeTextScalar} = "%.4g";
end

%% Get Values of er, ur, and thk
if isempty(thk)
    er = complex(O.erInitialValue, -O.erpInitialValue);
    ur = [];
    thk = O.thkInitialValue;
end

%% Check Inputs
if size(options.FormatString, 1) == 1
    options.FormatString = repmat(options.FormatString, length(thk), 1);
elseif size(options.FormatString, 1) ~= length(thk)
    error("'size(FormatString, 1)' must be either 1 or length(thk).");
end

if size(options.FormatString, 2) == 1
    options.FormatString = repmat(options.FormatString, 1, 5);
elseif size(options.FormatString, 2) ~= 5
    error("'size(FormatString, 2)' must be either 1 or 5.");
end

%% Add Limits for Each Parameter
if formatOptions.ShowLimits
    options.AdditionalText = strings(length(thk), 3, 2);
    
    for ii = 1:length(thk)
        % Lower Limit
        options.AdditionalText(ii, 1, 1) = sprintf(" [%s,", ...
            sprintf(options.FormatString(ii, 1), O.thkRange(1, ii)));
        
        options.AdditionalText(ii, 2, 1) = sprintf(" [%s - j%s,", ...
            sprintf(options.FormatString(ii, 2), O.erRange(1, ii)), ...
            sprintf(options.FormatString(ii, 3), O.erpRange(1, ii)));
        
        options.AdditionalText(ii, 3, 1) = sprintf(" [%s - j%s,", ...
            sprintf(options.FormatString(ii, 4), O.erRange(1, ii)), ...
            sprintf(options.FormatString(ii, 5), O.erpRange(1, ii)));
        
        % Upper Limit
        options.AdditionalText(ii, 1, 2) = sprintf("  %s]", ...
            sprintf(options.FormatString(ii, 1), O.thkRange(2, ii)));
        
        options.AdditionalText(ii, 2, 2) = sprintf("  %s - j%s]", ...
            sprintf(options.FormatString(ii, 2), O.erRange(2, ii)), ...
            sprintf(options.FormatString(ii, 3), O.erpRange(2, ii)));
        
        options.AdditionalText(ii, 3, 2) = sprintf("  %s - j%s]", ...
            sprintf(options.FormatString(ii, 4), O.erRange(2, ii)), ...
            sprintf(options.FormatString(ii, 5), O.erpRange(2, ii)));
    end
end

%% Customize Parameters that are Being Solved
if ~formatOptions.ShowInitialValues
    formatOptions.SolveParameterFormatString = sprintf(...
        formatOptions.SolveParameterFormatString, "%.0sx");
end

options.FormatString(O.thkLayersToSolve, 1) = compose(...
    formatOptions.SolveParameterFormatString, ...
    options.FormatString(O.thkLayersToSolve, 1));

options.FormatString(O.erLayersToSolve, 2) = compose(...
    formatOptions.SolveParameterFormatString, ...
    options.FormatString(O.erLayersToSolve, 2));

options.FormatString(O.erpLayersToSolve, 3) = compose(...
    formatOptions.SolveParameterFormatString, ...
    options.FormatString(O.erpLayersToSolve, 3));

%% Create String
optionsCell = namedargs2cell(options);
[varargout{1:nargout}] = nLayerForward.printStructure(er, ur, thk, ...
    optionsCell{1:end});

end
