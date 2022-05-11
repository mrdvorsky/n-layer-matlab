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
    options.ShowHyperlinks {mustBeNumericOrLogical} = true;
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

% Convert all "\n" to newlines.
structureString = compose(outputString);

%% Set MathText Formatted String
% figureString = replace(sprintf(outputString), "_ ", "  ");
figureString = sprintf(outputString);
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

% Remove bold from "additional text" lines.
numAdditionalLines = size(options.AdditionalText, 3);
indices = 5*(1:numAdditionalLines) + (1:numAdditionalLines).' - 1;
figureString(indices(:)) = strcat("\rm{", figureString(indices(:)), "}");
figureString(setdiff(1:length(figureString), indices(:))) = ...
    strcat("\bf{", figureString(setdiff(1:length(figureString), indices(:))), "}");

% Remove title line if empty.
if options.Title == ""
    figureString = figureString(2:end);
end

%% Add Hyperlink for Copy and Display
if options.ShowHyperlinks
    figureStringEscaped = replace(figureString, "\", "\\");
    hyperlink_base = "<a href=""matlab:%s"">%s</a>";
    
    commands_base = [...
        "hyp_fig = figure(Visible='off');", ...
        "hyp_ax = axes(hyp_fig, Position=[0, 0, 1, 1], XTick={}, YTick={}, ", ...
            "Color='none', XColor='none', YColor='none');", ...
        sprintf("hyp_string = sprintf(string('%s'));", strjoin(figureStringEscaped, "\n")), ...
        "hyp_text = text(hyp_ax, 0.5, 0.5, hyp_string, Units='normalized', ", ...
            "HorizontalAlignment='center', FontName='Monospaced', ", ...
            "FontWeight='bold', FontSize=10);", ...
        "hyp_fig.Position(3) = 1.01 * hyp_fig.Position(3) .* hyp_text.Extent(3);", ...
        "hyp_fig.Position(4) = 1.05 * hyp_fig.Position(4) .* hyp_text.Extent(4);", ...
        ];
    
    commands_show = [commands_base, "hyp_fig.Visible = 'on';"];
    commands_copy = [commands_base, "copygraphics(hyp_ax);", ...
        "close(hyp_fig);"];
    
    hyperlink_show = sprintf(hyperlink_base, strjoin(commands_show), "Show");
    hyperlink_copy = sprintf(hyperlink_base, strjoin(commands_copy), "Copy");
    
    % Change Title Line
    titleString = sprintf("%s (%s, %s)", ...
        options.Title, hyperlink_show, hyperlink_copy);
    
    structureStringArray = splitlines(structureString);
    structureStringArray(1) = pad(titleString, options.Width - 8 ...
        + strlength(hyperlink_show) + strlength(hyperlink_copy), "both");
    structureString = strjoin(structureStringArray, newline());
end

%% If Halfspace, Remove Backing Layer for Figure
if ~isfinite(thk(end))
    figureString(end) = "";
end

%% Print Output if Not Returning
if nargout == 0
    disp(strjoin(["", structureString, ""], newline()));
end

end

