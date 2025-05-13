function [varargout] = printStructureParameters(self, Parameters, Uncertainty, formatOptions, options)
%Display multilayer structure solver parameters.
% Prints the structure parameters for an nLayerInverse object. By default,
% uses the initial values from the nLayerInverse object. If a 'Parameters'
% struct is passed in (e.g., from the 'solveStructure' function), those
% values will be shown instead.
%
% Additionally, an 'Uncertainty' struct can be passed in to show
% uncertainty values (e.g., from the 'solveStructure' function).
%
% Example Usage:
%   NLsolver.printStructureParameters(Title="title", ShowLimits=true);
%   [Params, ~, Uncert] = NLsolver.solveStructure(...);
%   NLsolver.printStructureParams(Params);
%   NLsolver.printStructureParams(Params, Uncert);
%
%
% Author: Matt Dvorsky

arguments
    self nLayerInverse;

    Parameters(1, 1) {mustBeA(Parameters, "struct")} = struct();
    Uncertainty(1, 1) {mustBeA(Uncertainty, "struct")} = struct();

    formatOptions.ShowLimits(1, 1) logical = false;
    formatOptions.ShowInitialValues(1, 1) logical = true;
    formatOptions.SolveParameterFormatString {mustBeTextScalar} = "{%s}";

    options.Title {mustBeTextScalar} = "Solver Structure Parameters";
    options.BackingConductivity {mustBePositive} = inf;
    options.ThkUnitLabel {mustBeTextScalar} = "mm";
    options.ConductivityUnitLabel {mustBeTextScalar} = "S/mm";

    options.Width(1, 1) {mustBeInteger, mustBePositive} = 72;
    options.ThkWidth(1, 1) {mustBeInRange(options.ThkWidth, 0, 1)} = 0.3;
    options.ErWidth(1, 1) {mustBeInRange(options.ErWidth, 0, 1)} = 0.4;
    options.UrWidth(1, 1) {mustBeInRange(options.UrWidth, 0, 1)} = 0.3;
    options.IndentWidth(1, 1) {mustBeInteger, mustBeNonnegative} = 4;
    options.FormatString(:, :) {mustBeText} = "%.4g";
    options.ConductivityFormatString {mustBeTextScalar} = "%.4g";
end

%% Validate nLayerInverse Object
self.validate();

%% Get Values of er, ur, and thk
if isempty(fieldnames(Parameters))
    er = self.initialValue_er;
    ur = self.initialValue_ur;
    thk = self.initialValue_thk;
else
    er = Parameters.er;
    ur = Parameters.ur;
    thk = Parameters.thk;
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

%% Add Uncertainty Labels
options.AdditionalText = strings(numel(thk), 3, 0);
if ~isempty(fieldnames(Uncertainty))
    UncertaintyText = strings(length(thk), 3, 1);

    for ii = 1:length(thk)
        UncertaintyText(ii, 1, 1) = sprintf("     +-(%s)", ...
            sprintf(options.FormatString(ii, 1), Uncertainty.thk(ii)));

        UncertaintyText(ii, 2, 1) = sprintf("    +-(%s + j%s)", ...
            sprintf(options.FormatString(ii, 3), real(Uncertainty.er(ii))), ...
            sprintf(options.FormatString(ii, 4), imag(Uncertainty.er(ii))));

        UncertaintyText(ii, 3, 1) = sprintf("    +-(%s + j%s)", ...
            sprintf(options.FormatString(ii, 4), real(Uncertainty.ur(ii))), ...
            sprintf(options.FormatString(ii, 5), imag(Uncertainty.ur(ii))));
    end

    options.AdditionalText = cat(3, options.AdditionalText, UncertaintyText);
end

%% Add Limits for Each Parameter
if formatOptions.ShowLimits
    LimitsText = strings(length(thk), 3, 2);

    for ii = 1:length(thk)
        % Lower Limit
        LimitsText(ii, 1, 1) = sprintf(" %s,", ...
            sprintf(options.FormatString(ii, 1), self.rangeMin_thk(ii)));

        LimitsText(ii, 2, 1) = sprintf(" [%s - j%s,", ...
            sprintf(options.FormatString(ii, 2), self.rangeMin_erp(ii)), ...
            sprintf(options.FormatString(ii, 3), self.rangeMin_erpp(ii)));

        LimitsText(ii, 3, 1) = sprintf(" [%s - j%s,", ...
            sprintf(options.FormatString(ii, 4), self.rangeMin_urp(ii)), ...
            sprintf(options.FormatString(ii, 5), self.rangeMin_urpp(ii)));

        % Upper Limit
        LimitsText(ii, 1, 2) = sprintf("  %s]", ...
            sprintf(options.FormatString(ii, 1), self.rangeMax_thk(ii)));

        LimitsText(ii, 2, 2) = sprintf("  %s - j%s]", ...
            sprintf(options.FormatString(ii, 2), self.rangeMax_erp(ii)), ...
            sprintf(options.FormatString(ii, 3), self.rangeMax_erpp(ii)));

        LimitsText(ii, 3, 2) = sprintf("  %s - j%s]", ...
            sprintf(options.FormatString(ii, 4), self.rangeMax_urp(ii)), ...
            sprintf(options.FormatString(ii, 5), self.rangeMax_urpp(ii)));
    end

    options.AdditionalText = cat(3, options.AdditionalText, LimitsText);
end

%% Customize Parameters that are Being Solved
if ~formatOptions.ShowInitialValues
    formatOptions.SolveParameterFormatString = sprintf(...
        formatOptions.SolveParameterFormatString, "%.0sx");
end

options.FormatString(self.layersToSolve_thk, 1) = compose(...
    formatOptions.SolveParameterFormatString, ...
    options.FormatString(self.layersToSolve_thk, 1));

options.FormatString(self.layersToSolve_erp, 2) = compose(...
    formatOptions.SolveParameterFormatString, ...
    options.FormatString(self.layersToSolve_erp, 2));

options.FormatString(self.layersToSolve_erpp, 3) = compose(...
    formatOptions.SolveParameterFormatString, ...
    options.FormatString(self.layersToSolve_erpp, 3));

options.FormatString(self.layersToSolve_urp, 4) = compose(...
    formatOptions.SolveParameterFormatString, ...
    options.FormatString(self.layersToSolve_urp, 4));

options.FormatString(self.layersToSolve_urpp, 5) = compose(...
    formatOptions.SolveParameterFormatString, ...
    options.FormatString(self.layersToSolve_urpp, 5));

%% Create String
optionsCell = namedargs2cell(options);
[varargout{1:nargout}] = nLayer.printStructure(er, ur, thk, ...
    optionsCell{1:end});

end
