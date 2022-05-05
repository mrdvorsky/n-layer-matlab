function [structureString, figureString] = printStructure(er, ur, thk, options)
%PRINTSTRUCTURE Prints or returns a string showing nLayer structure.
%   Detailed explanation goes here
%
% Check sizes of inputs

arguments
    er(1, :);
    ur(1, :);
    thk(1, :) {mustBeNonempty};
    
    options.Title {mustBeTextScalar} = "Multilayer Structure Stackup";
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
    
    options.AdditionalText(:, 3, :) {mustBeText} = strings(1, 3, 0);
end

%% Check Inputs
[er, ur, thk] = nLayerForward.validateStructure(0, er, ur, thk, ...
    CheckStructureValues=false);

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

if size(options.AdditionalText, 1) == 1
    options.AdditionalText = repmat(options.AdditionalText, length(thk), 1);
elseif size(options.AdditionalText, 1) ~= length(thk)
    error("'size(AdditionalText, 1)' must be either 1 or length(thk).");
end

%% Create Helper Strings
titleString = pad(options.Title, options.Width, "both");
indentString = pad("", options.IndentWidth);
separatorString = strcat(pad("", options.Width, '_'), "\n");

if isfinite(thk(end))
    if isfinite(options.BackingConductivity)
        backingString = strcat(pad("", options.Width, '_'), "\n", ...
            pad(sprintf("  %s %s  ", ...
            sprintf(options.ConductivityFormatString, options.BackingConductivity), ...
            options.ConductivityUnitLabel), options.Width, "both", '/'));
    else
        backingString = strcat(pad("", options.Width, '_'), "\n", ...
            pad("", options.Width, '/'));
    end
else
    backingString = string(repmat('_ ', 1, round(0.5*options.Width)));
end

%% Calculate Table Widths
thkWidth = ceil(options.ThkWidth .* (options.Width - options.IndentWidth));
erWidth = ceil(options.ErWidth .* (options.Width - options.IndentWidth));
urWidth = ceil(options.UrWidth .* (options.Width - options.IndentWidth));

%% Create String for Each Layer
layerStrings = strings(length(thk), 1);
for ii = 1:length(thk)
    thkString = sprintf("thk%d = %s %s  ", ii, ...
        sprintf(options.FormatString(ii, 1), thk(ii)), ...
        options.ThkUnitLabel);
    
    erString = sprintf("er%d = %s - j%s  ", ii, ...
        sprintf(options.FormatString(ii, 2), real(er(ii))), ...
        sprintf(options.FormatString(ii, 3), abs(imag(er(ii)))));
    
    urString = sprintf("ur%d = %s - j%s", ii, ...
        sprintf(options.FormatString(ii, 4), real(ur(ii))), ...
        sprintf(options.FormatString(ii, 5), abs(imag(ur(ii)))));
    
    layerStrings(ii) = strcat(indentString, pad(thkString, thkWidth), ...
        pad(erString, erWidth), pad(urString, urWidth));
    
    
    % Add additional text if there is any
    for jj = 1:size(options.AdditionalText, 3)
        layerStrings(ii) = strcat(layerStrings(ii), "\n", indentString, ...
            pad(options.AdditionalText(ii, 1, jj), thkWidth), ...
            pad(options.AdditionalText(ii, 2, jj), erWidth), ...
            pad(options.AdditionalText(ii, 3, jj), urWidth));
    end
end

%% Combine Strings
outputString = strjoin([titleString, separatorString, ...
    strjoin(layerStrings, strcat("\n", separatorString, "\n")), ...
    backingString], "\n");

% Convert "\n" to newlines
structureString = sprintf(outputString);

%% Print Output if Not Returning
if nargout == 0
    fprintf(strcat("\n", outputString, "\n\n"));
end

%% Set MathText Formatted String
if nargout > 1
    figureString = replace(sprintf(outputString), "_ ", "  ");
    figureString = replace(figureString, "_", "\_");
    figureString = replace(figureString, "{", "\{");
    figureString = replace(figureString, "}", "\}");
    
    for ii = 1:length(thk)
        figureString = replace(figureString, ...
            sprintf("er%d", ii), sprintf("\\epsilon_{r{%d}}", ii));
        
        figureString = replace(figureString, ...
            sprintf("thk%d", ii), sprintf("thk_{%d}", ii));
        
        figureString = replace(figureString, ...
            sprintf("ur%d", ii), sprintf("\\mu_{r{%d}}", ii));
    end
    
    figureString = strtrim(splitlines(figureString));
end

end

