function [varargout] = nLayerViewer(er, ur, thk, NL, f, options)
%NLAYERVIEWER Function to open NLayer Viewer in a new figure window.
% This function allows users to compute multiple NLayer solvers and
% plot them to be viewed and compared. The sliders update in real-time
% allowing for a quick visual analysis of a material structure.
%
% Example Usage:
%   figure;
%   nlayerViewer(er, ur, thk, NL, f);
%   nlayerViewer(er, ur, thk, NL1, f1, NL2, f2, ...);
%   handles = nlayerViewer(er, ur, thk, NL1, f1, NL2, f2, ...);
%   handles = nlayerViewer(er, ur, thk, NL1, f1, ..., Legend=["1", ...]);
%
% Inputs:
%   er - Vector of default complex permittivities for each layer. These are
%       default values displayed for each layer.
%   ur - Vector of default complex permeabilities for each layer. These are
%       default values displayed for each layer.
%   thk - Vector of default layer thicknesses. Must have the same lenth as
%       er. Last value should be inf for infinite half-space.
%   NL - The nLayerForward solver object to use. This value is repeating
%       along with f.
%   f - Vector of frequencies. This value is repeating along with NL.
%
% Outputs:
%   handles - Optional output to return the handles containing data that is
%       stored within the UI figure.
%
% Named Arguments:
%   ShowLegend (true) - If true, legend will be enabled.
%   Legend - Array of strings that can be displayed along with the
%       various plots of different NLayer solvers. If not specified,
%       default names will be used.
%   ErBounds ([1, 10]) - Bounds of the dielectric constant (er). The
%       size of the variable is dependent on the number of layers
%       defined.
%   ErpBounds [0.01, 10]) - Bounds of the dielectric loss (erp). The size
%       of the variable is dependent on the number of layers defined.
%   UrBounds ([1, 10]) - Bounds of the magnetic constant (ur). The
%       size of the variable is dependent on the number of layers
%       defined.
%   UrpBounds [0.01, 10]) - Bounds of the magnetic loss (urp). The size
%       of the variable is dependent on the number of layers defined.
%   ThkBounds [0.1, 100]) - Bounds of the thickness (thk). The size of the
%       variable is dependent on the number of layers defined.
%   GamInterp (10) - The number of points the calculated reflection
%       coefficients is interpolated to.
%   PlotLineWidth (1.5) - The width of the plotted lines.
%   PlotDotWidth (1.5) - The width of the plotted dots.
%   FigureHandle (optional) - Handle of a figure to use. If not specified,
%       will use current active figure or create a new figure.
%   MainFigDim ([100, 100, 1000, 600]) - Main figure dimensions (pixels).
%       The array elements are [left, bottom, width, height].
%   PlotPanelSize (0.6) - Relative size of the plot panel to the main
%       figure window.
%   SliderXPos (0.15) - Position of sliders on the X-axis relative to
%       the size of the parent panel.
%   SliderYPos (0.15) - Position of sliders on the Y-axis relative to
%       the size of the parent panel.
%   SliderLength (0.6) - Length of the  sliders relative to the size of
%       the parent panel.
%   SliderWidth (0.12) - Width of the sliders relative to the size of
%       the parent panel.
%
% Author: Si Yuan Sim + Matt Dvorsky

arguments
    er  (1, :);
    ur  (1, :);
    thk (1, :) {mustBeNonempty};
end

arguments (Repeating)
    NL (1, 1) nLayerForward;
    f  (:, 1) {mustBeNonempty};
end

arguments
    % Figure options
    options.FigureHandle    (1, 1);
    options.MainFigDim      (1, 4) {mustBeReal} = [100, 100, 1000, 600];
    options.ShowLegend      (1, 1) {mustBeNumericOrLogical} = true;
    options.Legend          (:, 1) {mustBeText};
    options.PlotPanelSize   (1, 1) {mustBeReal} = 0.6;
    options.PanelFontSize   (1, 1) {mustBeReal} = 10;
    options.FigureColor = [0.94, 0.94, 0.94];

    % Specific Axis Sizes
    options.PlotAxisPosition   (1, 4) {mustBeReal} = [0.1, 0.03, 0.85, 0.9];
    options.StructureAxisSize  (1, 1) {mustBeReal};

    % Slider parameters
    options.SliderXPos      (1, 1) {mustBeReal} = 0.15;
    options.SliderYPos      (1, 1) {mustBeReal} = 0.15;
    options.SliderLength    (1, 1) {mustBeReal} = 0.6;
    options.SliderWidth     (1, 1) {mustBeReal} = 0.12;

    % Electrical property bounds
    options.ErBounds        (:, 2) {mustBePositive, mustBeFinite} = [1, 10];
    options.ErpBounds       (:, 2) {mustBeNonnegative, mustBeFinite} = [0.000, 1];
    options.UrBounds        (:, 2) {mustBePositive, mustBeFinite} = [1, 10];
    options.UrpBounds       (:, 2) {mustBeNonnegative, mustBeFinite} = [0.000, 1];
    options.ThkBounds       (:, 2) {mustBeNonnegative, mustBeFinite} = [0, 10];

    % NLayer settings
    options.NumFrequencyMarkers (1, 1) {mustBeInteger, mustBePositive} = 41;
    options.NumFrequencySamplesPerMarker (1, 1) {mustBeInteger, mustBePositive} = 20;

    % Plot settings
    options.PlotLineWidth   (1, 1) {mustBeReal} = 1;
    options.PlotMarkerSize  (1, 1) {mustBeReal} = 5;
    options.PlotMarkerType  {mustBeTextScalar} = ".";

    options.UpdateFunction = @(f, er, ur, thk) 0;
end

%% Check Inputs
if ~isfield(options, "StructureAxisSize")
    options.StructureAxisSize = 0.1 + 0.1*length(thk);
end

if (iscell(er))
    er = cell2mat(er);
end
if (iscell(ur))
    ur = cell2mat(ur);
end
if (iscell(thk))
    thk = cell2mat(thk);
end

if isempty(er)
    er = ones(size(thk));
elseif numel(er) == 1
    er = er + zeros(size(thk));
end

if isempty(ur)
    ur = ones(size(thk));
elseif numel(ur) == 1
    ur = ur + zeros(size(thk));
end

%% Create main figure and panels
if isfield(options, "FigureHandle")
    fig = options.FigureHandle;
else
    fig = gcf;
end
clf(fig);
fig.Position = options.MainFigDim;
fig.ToolBar = "none";

plotPanel = uipanel(fig, BackgroundColor=options.FigureColor, ...
    Position=[0, 0, options.PlotPanelSize, 1]);
sliderPanel = uipanel(fig, BackgroundColor=options.FigureColor, ...
    Position=[options.PlotPanelSize, 0, 1 - options.PlotPanelSize, 1]);

erPanel = uipanel(sliderPanel, Position=[0, 4/5, 1, 1/5], Tag="er", ...
    FontSize=options.PanelFontSize, Title="Dielectric Constant (er)");

erpPanel = uipanel(sliderPanel, Position=[0, 3/5, 1, 1/5], Tag="erp", ...
    FontSize=options.PanelFontSize, Title="Dielectric Loss (erp)");

urPanel = uipanel(sliderPanel, Position=[0, 2/5, 1, 1/5], Tag="ur", ...
    FontSize=options.PanelFontSize, Title="Magnetic Constant (urp)");

urpPanel = uipanel(sliderPanel, Position=[0, 1/5, 1, 1/5], Tag="urp", ...
    FontSize=options.PanelFontSize, Title="Magnetic Loss (urpp)");

thkPanel = uipanel(sliderPanel, Position=[0, 0, 1, 1/5], Tag="thk", ...
    FontSize=options.PanelFontSize, Title="Thickness (mm)");

handles.plotPanel = plotPanel;

%% Add Copy Figure Options to MenuBar
copyMenu = uimenu(fig, "Text", "Copy/Export Figure");

uimenu(copyMenu, "Text", "Copy Polar Plot", ...
    MenuSelectedFcn=@copyFigure);
uimenu(copyMenu, "Text", "Copy Structure Definition", ...
    MenuSelectedFcn=@copyStructure);
uimenu(copyMenu, "Text", "Copy Both", ...
    MenuSelectedFcn=@copyFigureAndStructure);

uimenu(copyMenu, "Text", "Export Polar Plot", ...
    MenuSelectedFcn=@exportFigure, Separator="on");
uimenu(copyMenu, "Text", "Export Structure Definition", ...
    MenuSelectedFcn=@exportStructure);
uimenu(copyMenu, "Text", "Export Both", ...
    MenuSelectedFcn=@exportFigureAndStructure);

%% Save initial material structure
numF = (options.NumFrequencyMarkers - 1) ...
    * options.NumFrequencySamplesPerMarker + 1;
for ii = 1:numel(f)
    f{ii} = linspace(min(f{ii}), max(f{ii}), numF);
end
handles.f = f;
handles.NL = cellfun(@copy, NL, UniformOutput=false);
handles.updateFunction = options.UpdateFunction;

numLayers = numel(thk);

%% Create Structure Description
structureAxis = axes(plotPanel, ...
    Position=[0, 0, 1, options.StructureAxisSize], ...
    TickLength=[0, 0], XTick={}, YTick={}, ...
    Color="none", XColor="none", YColor="none");

[~, structureString] = nLayer.printStructure(er, [], thk, Title="");

structureText = text(structureAxis, 0.5, 0.5, structureString, ...
    Units="normalized", HorizontalAlignment="center", ...
    FontName="Monospaced", FontWeight="bold", FontSize=9);

handles.structureText = structureText;
handles.structureAxis = structureAxis;

%% Create Polar Plot
plotAxisPosition = options.PlotAxisPosition;
plotAxisPosition(2) = 1 - (1 - plotAxisPosition(2)) ...
    .* (1 - options.StructureAxisSize);
plotAxisPosition(4) = plotAxisPosition(4) ...
    .* (1 - options.StructureAxisSize);
plotAxis = axes(plotPanel, Position=plotAxisPosition, ...
    ButtonDownFcn=@buttonClickFunction);

[h1, h2, h3] = zplane([], [], plotAxis);
h1.HandleVisibility = "off";
h2.HandleVisibility = "off";
h3.HandleVisibility = "off";
xlabel(plotAxis, "");
ylabel(plotAxis, "");
title(plotAxis, "");

xlim(plotAxis, [-1.1, 1.1]);
ylim(plotAxis, [-1.1, 1.1]);
hold(plotAxis, "on");

gamPlot = cell(numel(NL), 1);
gamFitPlot = cell(numel(NL), 1);

% Obtain initial parameters and calculate initial values
for ii = 1:numel(NL)
    gam = NL{ii}.calculate(f{ii}, er, ur, thk);

    plotLabels = NL{ii}.getOutputLabels();
    gamFitPlot{ii} = cell(size(gam(:, :), 2), 1);

    for pp = 1:size(gam(:, :), 2)
        gamFitPlot{ii}{pp} = plot(plotAxis, gam(:, pp), ...
            Linewidth=options.PlotLineWidth, ...
            DisplayName=plotLabels(pp), ...
            Marker=options.PlotMarkerType, ...
            MarkerSize=options.PlotMarkerSize, ...
            MarkerIndices=(0:options.NumFrequencyMarkers - 1) * options.NumFrequencySamplesPerMarker + 1);
    end
end

handles.markerPlot = plot(-1, -1, ".k", MarkerSize=8, ...
    HitTest="off", DisplayName=sprintf("f = %g GHz", f{1}(1)));
handles.markerFreq = f{1}(1);

for ii = 1:numel(plotAxis.Children)
    plotAxis.Children(ii).HitTest = "off";
end

if options.ShowLegend
    if isfield(options, "Legend")
        legend(plotAxis, options.Legend);
    else
        legend(plotAxis);
    end
end

hold(plotAxis, "off");

handles.plotAxis = plotAxis;
handles.gamPlot = gamPlot;
handles.gamFitPlot = gamFitPlot;

%% Figure Size Change Callback
fig.SizeChangedFcn = @figureResizeCallback;

%% Create sliders
createSlider = @(panel, ind, tag) uicontrol(Parent=panel, ...
    Style="slider", Units="normalized", ...
    Position=[options.SliderXPos, 0.75 - options.SliderYPos*(ind - 1), ...
    options.SliderLength, options.SliderWidth], Tag=tag);

uiSliders.erSliders =  cell(numLayers, 1);
uiSliders.erpSliders = cell(numLayers, 1);
uiSliders.urSliders =  cell(numLayers, 1);
uiSliders.urpSliders = cell(numLayers, 1);
uiSliders.thkSliders = cell(numLayers, 1);
uiSliders.erRange =  options.ErBounds  + zeros(numLayers, 2);
uiSliders.erpRange = options.ErpBounds + zeros(numLayers, 2);
uiSliders.urRange =  options.UrBounds  + zeros(numLayers, 2);
uiSliders.urpRange = options.UrpBounds + zeros(numLayers, 2);
uiSliders.thkRange = options.ThkBounds + zeros(numLayers, 2);

for ii = 1:numLayers
    % Create sliders
    uiSliders.erSliders{ii}  = createSlider(erPanel,  ii, sprintf("er-%d", ii));
    uiSliders.erpSliders{ii} = createSlider(erpPanel, ii, sprintf("erp-%d", ii));
    uiSliders.urSliders{ii}  = createSlider(urPanel,  ii, sprintf("ur-%d", ii));
    uiSliders.urpSliders{ii} = createSlider(urpPanel, ii, sprintf("urp-%d", ii));
    uiSliders.thkSliders{ii} = createSlider(thkPanel, ii, sprintf("thk-%d", ii));

    % Preset initial values
    erSliderLoc = in_lerp(real(er(ii)), uiSliders.erRange(ii, :));
    if erSliderLoc > 1
        uiSliders.erSliders{ii}.Value = 1;
    elseif erSliderLoc < 0
        uiSliders.erSliders{ii}.Value = 0;
    else
        uiSliders.erSliders{ii}.Value = erSliderLoc;
    end

    erpSliderLoc = in_lerp((-imag(er(ii))), uiSliders.erpRange(ii, :));
    if erpSliderLoc > 1
        uiSliders.erpSliders{ii}.Value = 1;
    elseif erpSliderLoc < 0
        uiSliders.erpSliders{ii}.Value = 0;
    else
        uiSliders.erpSliders{ii}.Value = erpSliderLoc;
    end

    urSliderLoc = in_lerp(real(ur(ii)), uiSliders.urRange(ii, :));
    if urSliderLoc > 1
        uiSliders.urSliders{ii}.Value = 1;
    elseif urSliderLoc < 0
        uiSliders.urSliders{ii}.Value = 0;
    else
        uiSliders.urSliders{ii}.Value = urSliderLoc;
    end

    urpSliderLoc = in_lerp((-imag(ur(ii))), uiSliders.urpRange(ii, :));
    if urpSliderLoc > 1
        uiSliders.urpSliders{ii}.Value = 1;
    elseif urpSliderLoc < 0
        uiSliders.urpSliders{ii}.Value = 0;
    else
        uiSliders.urpSliders{ii}.Value = urpSliderLoc;
    end

    thkSliderLoc = in_lerp(thk(ii), uiSliders.thkRange(ii, :));
    if thkSliderLoc > 1
        uiSliders.thkSliders{ii}.Value = 1;
    elseif thkSliderLoc < 0
        uiSliders.thkSliders{ii}.Value = 0;
    else
        uiSliders.thkSliders{ii}.Value = thkSliderLoc;
    end
end

% Create checkbox for infinite half plane in thk panel
handles.isInfHalfPlane = uicontrol(Style="checkbox", Parent=thkPanel, ...
    String="Infinite Half-Plane (Bottom Layer)", Units="normalized", ...
    Position=[0, 0.75 - 0.15*numLayers, 0.5, 0.1], ...
    FontSize=8, HorizontalAlignment="right", CallBack=@halfPlaneValueChange);

%% Create slider labels
createTopLabel = @(panel) uicontrol(Style="text", String="Top Layer", ...
    Parent=panel, Units="normalized", Position=[0, 0.9, 1, 0.1], ...
    HorizontalAlignment="left");

createBottomLabel = @(panel) uicontrol(Style="text", String="Bottom Layer", ...
    Parent=panel, Units="normalized", ...
    Position=[0, 0.75 - options.SliderYPos*numLayers, 0.2, 0.1], ...
    HorizontalAlignment="left");

createTopLabel(erPanel);
createBottomLabel(erPanel);
createTopLabel(erpPanel);
createBottomLabel(erpPanel);
createTopLabel(urPanel);
createBottomLabel(urPanel);
createTopLabel(urpPanel);
createBottomLabel(urpPanel);
createTopLabel(thkPanel);

%% Create edit field
uiEditField.erLB = cell(numLayers, 1);
uiEditField.erUB = cell(numLayers, 1);
uiEditField.erCV = cell(numLayers, 1);
uiEditField.erpLB = cell(numLayers, 1);
uiEditField.erpUB = cell(numLayers, 1);
uiEditField.erpCV = cell(numLayers, 1);
uiEditField.urLB = cell(numLayers, 1);
uiEditField.urUB = cell(numLayers, 1);
uiEditField.urCV = cell(numLayers, 1);
uiEditField.urpLB = cell(numLayers, 1);
uiEditField.urpUB = cell(numLayers, 1);
uiEditField.urpCV = cell(numLayers, 1);
uiEditField.thkLB = cell(numLayers, 1);
uiEditField.thkUB = cell(numLayers, 1);
uiEditField.thkCV = cell(numLayers, 1);

% Create edit fields
for ii = 1:numLayers
    uiEditField.erLB{ii} = LBValueField(erPanel, ii, lerp(0, uiSliders.erRange(ii, :)));
    uiEditField.erUB{ii} = UBValueField(erPanel, ii, lerp(1, uiSliders.erRange(ii, :)));
    uiEditField.erCV{ii} = CVValueField(erPanel, ii, real(er(ii)));

    uiEditField.erpLB{ii} = LBValueField(erpPanel, ii, lerp(0, uiSliders.erpRange(ii, :)));
    uiEditField.erpUB{ii} = UBValueField(erpPanel, ii, lerp(1, uiSliders.erpRange(ii, :)));
    uiEditField.erpCV{ii} = CVValueField(erpPanel, ii, (-imag(er(ii))));

    uiEditField.urLB{ii} = LBValueField(urPanel, ii, lerp(0, uiSliders.urRange(ii, :)));
    uiEditField.urUB{ii} = UBValueField(urPanel, ii, lerp(1, uiSliders.urRange(ii, :)));
    uiEditField.urCV{ii} = CVValueField(urPanel, ii, real(ur(ii)));

    uiEditField.urpLB{ii} = LBValueField(urpPanel, ii, lerp(0, uiSliders.urpRange(ii, :)));
    uiEditField.urpUB{ii} = UBValueField(urpPanel, ii, lerp(1, uiSliders.urpRange(ii, :)));
    uiEditField.urpCV{ii} = CVValueField(urpPanel, ii, (-imag(ur(ii))));

    uiEditField.thkLB{ii} = LBValueField(thkPanel, ii, lerp(0, uiSliders.thkRange(ii, :)));
    uiEditField.thkUB{ii} = UBValueField(thkPanel, ii, lerp(1, uiSliders.thkRange(ii, :)));
    uiEditField.thkCV{ii} = CVValueField(thkPanel, ii, thk(ii));
end

handles.uiEditField = uiEditField;

%% Create slider function handle wrapper
valueChange = @(hObject, eventdata) sliderValueChanged(hObject, eventdata);

for ii = 1:numLayers
    % Set callbacks
    uiSliders.erSliders{ii}.Callback = valueChange;
    uiSliders.erpSliders{ii}.Callback = valueChange;
    uiSliders.urSliders{ii}.Callback = valueChange;
    uiSliders.urpSliders{ii}.Callback = valueChange;
    uiSliders.thkSliders{ii}.Callback = valueChange;

    % Set listeners
    addlistener(uiSliders.erSliders{ii}, "Value", "PostSet", ...
        @(~, eventdata) sliderValueChanged(eventdata.AffectedObject, eventdata));
    addlistener(uiSliders.erpSliders{ii}, "Value", "PostSet", ...
        @(~, eventdata) sliderValueChanged(eventdata.AffectedObject, eventdata));
    addlistener(uiSliders.urSliders{ii}, "Value", "PostSet", ...
        @(~, eventdata) sliderValueChanged(eventdata.AffectedObject, eventdata));
    addlistener(uiSliders.urpSliders{ii}, "Value", "PostSet", ...
        @(~, eventdata) sliderValueChanged(eventdata.AffectedObject, eventdata));
    addlistener(uiSliders.thkSliders{ii}, "Value", "PostSet", ...
        @(~, eventdata) sliderValueChanged(eventdata.AffectedObject, eventdata));
end

handles.uiSliders = uiSliders;
handles.options = options;

% Handles is an optional output parameter
if nargout == 1
    varargout{1} = handles;
end

% Store data into figure
guidata(fig, handles);

% Trigger callback to reset display
CVFieldChanged(uiEditField.thkCV{1}, []);

end

%% Slider value change function
function sliderValueChanged(hObject, ~)
handles = guidata(hObject);

% Extract value from slider and current value edit field
[er, erp, ur, urp, thk] = valueReader(handles);
[er_slider, erp_slider, ur_slider, urp_slider, thk_slider] = ...
    valueExtractor(handles);

panel = extractBefore(hObject.Tag, "-");
layer = str2double(extractAfter(hObject.Tag, "-"));

% Update the value being changed in the appropriate edit field
switch panel
    case "er"
        handles.uiEditField.erCV{layer}.String = er_slider(layer);
        er(layer) = er_slider(layer);
    case "erp"
        handles.uiEditField.erpCV{layer}.String = erp_slider(layer);
        erp(layer) = erp_slider(layer);
    case "ur"
        handles.uiEditField.urCV{layer}.String = ur_slider(layer);
        ur(layer) = ur_slider(layer);
    case "urp"
        handles.uiEditField.urpCV{layer}.String = urp_slider(layer);
        urp(layer) = urp_slider(layer);
    case "thk"
        handles.uiEditField.thkCV{layer}.String = thk_slider(layer);
        thk(layer) = thk_slider(layer);
end

handles = plotGam(handles);

% Store data into figure
guidata(hObject, handles);

% Force redraw axes
drawnow;

end

%% Lower bound edit field creation
function LBField = LBValueField(panel, ind, initVal)
LBField = uicontrol(Style="edit", Parent=panel, Units="Normalized", ...
    CallBack=@LBFieldChanged, Tag=num2str(ind), String=num2str(initVal), ...
    Position=[0.05, 0.75 - 0.15*(ind - 1), 0.1, 0.12]);

% uicontrol(LBField);
end

%% Lower bound edit field callback
function LBFieldChanged(hObject, ~)
handles = guidata(hObject);

panel = hObject.Parent.Tag;
layer = str2double(hObject.Tag);
lowerBound = str2double(hObject.String);

if ~isnan(lowerBound)
    switch panel
        case "er"
            currentValue = str2double(handles.uiEditField.erCV{layer}.String);
            if lowerBound < 0
                hObject.String = num2str(1);
                handles.uiSliders.erRange(layer,1) = 0;
            elseif lowerBound <= currentValue
                handles.uiSliders.erRange(layer,1) = lowerBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.erRange(layer,1) = currentValue;
            end
            guidata(hObject, handles);
        case "erp"
            currentValue = str2double(handles.uiEditField.erpCV{layer}.String);
            if lowerBound <= currentValue
                handles.uiSliders.erpRange(layer,1) = lowerBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.erpRange(layer,1) = currentValue;
            end
            guidata(hObject, handles);
        case "ur"
            currentValue = str2double(handles.uiEditField.urCV{layer}.String);
            if lowerBound < 0
                hObject.String = num2str(1);
                handles.uiSliders.urRange(layer,1) = 0;
            elseif lowerBound <= currentValue
                handles.uiSliders.urRange(layer,1) = lowerBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.urRange(layer,1) = currentValue;
            end
            guidata(hObject, handles);
        case "urp"
            currentValue = str2double(handles.uiEditField.urpCV{layer}.String);
            if lowerBound <= currentValue
                handles.uiSliders.urpRange(layer,1) = lowerBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.urpRange(layer,1) = currentValue;
            end
            guidata(hObject, handles);
        case "thk"
            currentValue = str2double(handles.uiEditField.thkCV{layer}.String);
            if lowerBound <= currentValue
                handles.uiSliders.thkRange(layer,1) = lowerBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.thkRange(layer,1) = currentValue;
            end
            guidata(hObject, handles);
    end
else
    switch panel
        case "er"
            hObject.String = num2str(handles.uiSliders.erRange(layer,1));
        case "erp"
            hObject.String = num2str(handles.uiSliders.erpRange(layer,1));
        case "ur"
            hObject.String = num2str(handles.uiSliders.urRange(layer,1));
        case "urp"
            hObject.String = num2str(handles.uiSliders.urpRange(layer,1));
        case "thk"
            hObject.String = num2str(handles.uiSliders.thkRange(layer,1));
    end
end

% Store data into figure
guidata(hObject, handles);

end

%% Upper bound edit field and callback
function limitField = UBValueField(panel, ind, initVal)
limitField = uicontrol(Style="edit", Parent=panel, Units="Normalized", ...
    CallBack=@UBFieldChanged, Tag=num2str(ind), String=num2str(initVal), ...
    Position=[0.745, 0.75 - 0.15*(ind - 1), 0.1, 0.12]);

uicontrol(limitField);
end

function UBFieldChanged(hObject, ~)
handles = guidata(hObject);

panel = hObject.Parent.Tag;
layer = str2double(hObject.Tag);
upperBound = str2double(hObject.String);

if ~isnan(upperBound)
    switch panel
        case "er"
            currentValue = str2double(handles.uiEditField.erCV{layer}.String);

            if upperBound >= currentValue
                handles.uiSliders.erRange(layer,2) = upperBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.erRange(layer,2) = currentValue;
            end

            guidata(hObject, handles);
            handles.uiSliders.erSliders{layer}.Value = in_lerp(currentValue, handles.uiSliders.erRange(layer, :));
        case "erp"
            currentValue = str2double(handles.uiEditField.erpCV{layer}.String);

            if upperBound >= currentValue
                handles.uiSliders.erpRange(layer,2) = upperBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.erpRange(layer,2) = currentValue;
            end

            guidata(hObject, handles);
            handles.uiSliders.erpSliders{layer}.Value = in_lerp(currentValue, handles.uiSliders.erpRange(layer, :));
        case "ur"
            currentValue = str2double(handles.uiEditField.urCV{layer}.String);

            if upperBound >= currentValue
                handles.uiSliders.urRange(layer,2) = upperBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.urRange(layer,2) = currentValue;
            end

            guidata(hObject, handles);
            handles.uiSliders.urSliders{layer}.Value = in_lerp(currentValue, handles.uiSliders.urRange(layer, :));
        case "urp"
            currentValue = str2double(handles.uiEditField.urpCV{layer}.String);

            if upperBound >= currentValue
                handles.uiSliders.urpRange(layer,2) = upperBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.urpRange(layer,2) = currentValue;
            end

            guidata(hObject, handles);
            handles.uiSliders.urpSliders{layer}.Value = in_lerp(currentValue, handles.uiSliders.urpRange(layer, :));
        case "thk"
            currentValue = str2double(handles.uiEditField.thkCV{layer}.String);

            if upperBound >= currentValue
                handles.uiSliders.thkRange(layer,2) = upperBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.thkRange(layer,2) = currentValue;
            end

            guidata(hObject, handles);
            handles.uiSliders.thkSliders{layer}.Value = in_lerp(currentValue, handles.uiSliders.thkRange(layer, :));
    end
else
    switch panel
        case "er"
            hObject.String = num2str(handles.uiSliders.erRange(layer,2));
        case "erp"
            hObject.String = num2str(handles.uiSliders.erpRange(layer,2));
        case "ur"
            hObject.String = num2str(handles.uiSliders.urRange(layer,2));
        case "urp"
            hObject.String = num2str(handles.uiSliders.urpRange(layer,2));
        case "thk"
            hObject.String = num2str(handles.uiSliders.thkRange(layer,2));
    end
end

% Store data into figure
guidata(hObject, handles);

end

%% Current value edit field creation
function CVField = CVValueField(panel, ind, initVal)
CVField = uicontrol(Style="edit", Parent=panel, Units="Normalized", ...
    CallBack=@CVFieldChanged, Tag=num2str(ind), String=num2str(initVal), ...
    Position=[0.87, 0.75 - 0.15*(ind - 1), 0.1, 0.12]);

uicontrol(CVField);
end

%% Current value edit field callback
function CVFieldChanged(hObject, ~)
handles = guidata(hObject);

uiSliders = handles.uiSliders;

panel = hObject.Parent.Tag;
layer = str2double(hObject.Tag);
currentValue = str2double(hObject.String);
isOutsideBounds = 0;

if ~isnan(currentValue) || isreal(currentValue)
    switch panel
        case "er"
            if currentValue < 1
                hObject.String = num2str(lerp(uiSliders.erSliders{layer}.Value, uiSliders.erRange(layer, :)));
            end

            sliderValue = in_lerp(currentValue, uiSliders.erRange(layer, :));

            if sliderValue >= 0 && sliderValue <= 1
                uiSliders.erSliders{layer}.Value = sliderValue;
            elseif sliderValue < 0
                uiSliders.erSliders{layer}.Value = 0;
                isOutsideBounds = 1;
            elseif sliderValue > 1
                uiSliders.erSliders{layer}.Value = 1;
                isOutsideBounds = 1;
            else
                hObject.String = num2str(lerp(uiSliders.erSliders{layer}.Value, uiSliders.erRange(layer, :)));
            end
        case "erp"
            sliderValue = in_lerp(currentValue, uiSliders.erpRange(layer, :));

            if sliderValue >= 0 && sliderValue <= 1
                uiSliders.erpSliders{layer}.Value = sliderValue;
            elseif sliderValue < 0
                uiSliders.erpSliders{layer}.Value = 0;
                isOutsideBounds = 1;
            elseif sliderValue > 1
                uiSliders.erpSliders{layer}.Value = 1;
                isOutsideBounds = 1;
            else
                hObject.String = num2str(lerp(uiSliders.erpSliders{layer}.Value, uiSliders.erpRange(layer, :)));
            end
        case "ur"
            if currentValue < 1
                hObject.String = num2str(lerp(uiSliders.urSliders{layer}.Value, uiSliders.urRange(layer, :)));
            end

            sliderValue = in_lerp(currentValue, uiSliders.urRange(layer, :));

            if sliderValue >= 0 && sliderValue <= 1
                uiSliders.urSliders{layer}.Value = sliderValue;
            elseif sliderValue < 0
                uiSliders.urSliders{layer}.Value = 0;
                isOutsideBounds = 1;
            elseif sliderValue > 1
                uiSliders.urSliders{layer}.Value = 1;
                isOutsideBounds = 1;
            else
                hObject.String = num2str(lerp(uiSliders.urSliders{layer}.Value, uiSliders.urRange(layer, :)));
            end
        case "urp"
            sliderValue = in_lerp(currentValue, uiSliders.urpRange(layer, :));

            if sliderValue >= 0 && sliderValue <= 1
                uiSliders.urpSliders{layer}.Value = sliderValue;
            elseif sliderValue < 0
                uiSliders.urpSliders{layer}.Value = 0;
                isOutsideBounds = 1;
            elseif sliderValue > 1
                uiSliders.urpSliders{layer}.Value = 1;
                isOutsideBounds = 1;
            else
                hObject.String = num2str(lerp(uiSliders.urpSliders{layer}.Value, uiSliders.urpRange(layer, :)));
            end
        case "thk"
            sliderValue = in_lerp(currentValue, uiSliders.thkRange(layer, :));

            if sliderValue >= 0 && sliderValue <= 1
                uiSliders.thkSliders{layer}.Value = sliderValue;
            elseif sliderValue < 0
                uiSliders.thkSliders{layer}.Value = 0;
                isOutsideBounds = 1;
            elseif sliderValue > 1
                uiSliders.thkSliders{layer}.Value = 1;
                isOutsideBounds = 1;
            else
                hObject.String = num2str(lerp(uiSliders.thkSliders{layer}.Value, uiSliders.thkRange(layer, :)));
            end
    end
else
    switch panel
        case "er"
            hObject.String = num2str(lerp(uiSliders.erSliders{layer}.Value, uiSliders.erRange(layer, :)));
        case "erp"
            hObject.String = num2str(lerp(uiSliders.erpSliders{layer}.Value, uiSliders.erpRange(layer, :)));
        case "ur"
            hObject.String = num2str(lerp(uiSliders.urSliders{layer}.Value, uiSliders.urRange(layer, :)));
        case "urp"
            hObject.String = num2str(lerp(uiSliders.urpSliders{layer}.Value, uiSliders.urpRange(layer, :)));
        case "thk"
            hObject.String = num2str(lerp(uiSliders.thkSliders{layer}.Value, uiSliders.thkRange(layer, :)));
    end
end

if isOutsideBounds
    hObject.String = currentValue;
end

handles = plotGam(handles);

handles.uiSliders = uiSliders;

% Store data into figure
guidata(hObject, handles);

end

%% Infinite half plane checkbox value change callback
function halfPlaneValueChange(hObject, ~)
handles = guidata(hObject);

isInfHalfPlane = handles.isInfHalfPlane.Value;
if isInfHalfPlane
    handles.uiSliders.thkSliders{end}.Enable = "off";
    handles.uiEditField.thkCV{end}.Enable = "off";
    handles.uiEditField.thkLB{end}.Enable = "off";
    handles.uiEditField.thkUB{end}.Enable = "off";
else
    handles.uiSliders.thkSliders{end}.Enable = "on";
    handles.uiEditField.thkCV{end}.Enable = "on";
    handles.uiEditField.thkLB{end}.Enable = "on";
    handles.uiEditField.thkUB{end}.Enable = "on";
end

handles = plotGam(handles);
guidata(hObject, handles);

end

%% Interpolate and inverse interpolation of values
function [v] = lerp(m, range)
v = range(1) + m.*(range(end) - range(1));
end

function [m] = in_lerp(v, range)
m = (v - range(1))./(range(end) - range(1));
end

function [v_arr] = lerp_arr(m, range)
v_arr = range(:, 1) + m.*(range(:, end) - range(:, 1));
end

%% Value extractor
function [er, erp, ur, urp, thk] = valueExtractor(handles)
valueExtracted = @(s) s.Value;
er =  lerp_arr(cellfun(valueExtracted, handles.uiSliders.erSliders), handles.uiSliders.erRange);
erp = lerp_arr(cellfun(valueExtracted, handles.uiSliders.erpSliders), handles.uiSliders.erpRange);
ur =  lerp_arr(cellfun(valueExtracted, handles.uiSliders.urSliders), handles.uiSliders.urRange);
urp = lerp_arr(cellfun(valueExtracted, handles.uiSliders.urpSliders), handles.uiSliders.urpRange);
thk = lerp_arr(cellfun(valueExtracted, handles.uiSliders.thkSliders), handles.uiSliders.thkRange);
end

%% Value reader from current value edit field
function [erp, erpp, urp, urpp, thk] = valueReader(handles)
% Read current value from edit field which is string stored in cell array
valueRead = @(s) s.String;
er_str = cellfun(valueRead, handles.uiEditField.erCV, UniformOutput=false);
erp_str = cellfun(valueRead, handles.uiEditField.erpCV, UniformOutput=false);
ur_str = cellfun(valueRead, handles.uiEditField.urCV, UniformOutput=false);
urp_str = cellfun(valueRead, handles.uiEditField.urpCV, UniformOutput=false);
thk_str = cellfun(valueRead, handles.uiEditField.thkCV, UniformOutput=false);

% Convert string to number
erp = cellfun(@str2double, er_str);
erpp = cellfun(@str2double, erp_str);
urp = cellfun(@str2double, ur_str);
urpp = cellfun(@str2double, urp_str);
thk = cellfun(@str2double, thk_str);

if handles.isInfHalfPlane.Value
    thk(end) = inf;
end

end

function [er, ur, thk] = valueReaderComplex(handles)
[erp, erpp, urp, urpp, thk] = valueReader(handles);
er = num2cell(erp - 1j*erpp);
ur = num2cell(urp - 1j*urpp);
thk = num2cell(thk);
end

%% Figure Resize Callback
function figureResizeCallback(fig, ~)
handles = guidata(fig);
figMinRatio = min(fig.Position([3, 4]) ./ [1000, 600]);
handles.structureText.FontSize = round(9 .* figMinRatio);
end

%% Copy Figure Callbacks
function copyFigure(fig, ~)
handles = guidata(fig.Parent.Parent);
copygraphics(handles.plotAxis);
end

function copyStructure(fig, ~)
handles = guidata(fig.Parent.Parent);
copygraphics(handles.structureAxis);
end

function copyFigureAndStructure(fig, ~)
handles = guidata(fig.Parent.Parent);
copygraphics(handles.plotPanel);
end

%% Export Figure Callbacks
function exportFigure(fig, ~)
handles = guidata(fig.Parent.Parent);
newFig = figure;
newPlotAxis = copyobj(handles.plotAxis, newFig);
newPlotAxis.Position = [0.13, 0.11, 0.775, 0.815];

if isprop(handles.plotAxis, "Legend")
    legend(newPlotAxis);
end
end

function exportStructure(fig, ~)
handles = guidata(fig.Parent.Parent);
newFig = figure();
newStructureAxis = copyobj(handles.structureAxis, newFig);
newStructureAxis.Children(1).FontSize = 9;
newStructureAxis.Position = [0.1, 0.1, 0.8, 0.8];
end

function exportFigureAndStructure(fig, ~)
handles = guidata(fig.Parent.Parent);
newFig = figure(Position=[680, 200, 560, 600]);
newPlotAxis = copyobj(handles.plotAxis, newFig);
newStructureAxis = copyobj(handles.structureAxis, newFig);
newStructureAxis.Children(1).FontSize = 9;

if isprop(handles.plotAxis, "Legend")
    legend(newPlotAxis);
end
end

%% Calculate and plot S-Parameters
function handles = plotGam(handles)
[er, ur, thk] = valueReaderComplex(handles);

for ii = 1:size(handles.NL, 2)
    gam = handles.NL{ii}.calculate(handles.f{ii}, er, ur, thk);

    for pp = 1:size(gam(:, :), 2)
        handles.gamFitPlot{ii}{pp}.XData = real(gam(:, pp));
        handles.gamFitPlot{ii}{pp}.YData = imag(gam(:, pp));
    end
end

gamF = handles.NL{1}.calculate(handles.markerFreq, er, ur, thk);
handles.markerPlot.XData = real(gamF);
handles.markerPlot.YData = imag(gamF);
handles.updateFunction(handles.markerFreq, er, ur, thk);

[~, handles.structureText.String] = ...
    nLayer.printStructure(er, ur, thk, Title="");

end

%% Button Click CallBack
function [pIndFrac] = nearest_point(line_x, line_y, xp, yp)
    ux = line_x(1:end-1) - xp;
    uy = line_y(1:end-1) - yp;
    vx = line_x(2:end) - line_x(1:end-1);
    vy = line_y(2:end) - line_y(1:end-1);

    t = max(0, min(1, -(ux.*vx + uy.*vy) ./ (vx.*vx + vy.*vy)));
    cx = (1 - t).*line_x(1:end-1) + t.*line_x(2:end);
    cy = (1 - t).*line_y(1:end-1) + t.*line_y(2:end);

    [~, pInd] = min(hypot(cx - xp, cy - yp));
    pIndFrac = pInd + t(pInd);

    % [~, pInd] = min(hypot(line_x - xp, line_y - yp));
end

function buttonClickFunction(hObject, ~)
    handles = guidata(hObject);
    xy = hObject.CurrentPoint(1, 1:2);

    line = handles.gamFitPlot{1}{1};
    pInd = nearest_point(line.XData, line.YData, xy(1), xy(2));
    f = interp1(handles.f{1}, pInd);
    
    [er, ur, thk] = valueReaderComplex(handles);
    gamF = handles.NL{1}.calculate(f, er, ur, thk);

    handles.markerPlot.XData = real(gamF);
    handles.markerPlot.YData = imag(gamF);
    handles.markerFreq = f;
    handles.markerPlot.DisplayName = sprintf("f = %g GHz", f);
    handles.updateFunction(handles.markerFreq, er, ur, thk);

    guidata(hObject, handles);
end



