function [varargout] = nLayerViewer(er, thk, NL, f, options)
%NLAYERVIEWER Function to open NLayer Viewer in a new figure window.
% This function allows users to compute multiple NLayer solvers and
% plot them to be viewed and compared. The sliders update in real-time
% allowing for a quick visual analysis of a material structure.
%
% Example Usage:
%   figure;
%   nlayerViewer(er, thk, NL, f);
%   nlayerViewer(er, thk, NL1, f1, NL2, f2, ...);
%   handles = nlayerViewer(er, thk, NL1, f1, NL2, f2, ...);
%   handles = nlayerViewer(er, thk, NL1, f1, ..., Legend=["1", ...]);
%
% Inputs:
%   er - Vector of default complex permittivities for each layer. These are
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
%   Legend - Array of strings that can be displayed along with the
%       various plots of different NLayer solvers.
%   ErBounds ([1, 10]) - Bounds of the dielectric constant (er). The
%       size of the variable is dependent on the number of layers
%       defined.
%   ErpBounds [-2, 1]) - Bounds of the dielectric loss (erp). The value
%       is defined as the integer power to 10. The size of the variable
%       is dependent on the number of layers defined.
%   ThkBounds [-1, 2]) - Bounds of the thickness (thk). The value
%       is defined as the integer power to 10. The size of the variable
%       is dependent on the number of layers defined.
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
    thk (1, :);
end

arguments (Repeating)
    NL (1, 1);
    f  (:, 1);
end

arguments
    % Figure options
    options.FigureHandle    (1, 1);
    options.MainFigDim      (1, 4) {mustBeNumeric} = [100, 100, 1000, 600];
    options.Legend          (:, 1) {mustBeText};
    options.PlotPanelSize   (1, 1) {mustBeNumeric} = 0.6;
    options.PanelFontSize   (1, 1) {mustBeNumeric} = 10;
    
    % Slider parameters
    options.SliderXPos      (1, 1) {mustBeNumeric} = 0.15;
    options.SliderYPos      (1, 1) {mustBeNumeric} = 0.15;
    options.SliderLength    (1, 1) {mustBeNumeric} = 0.6;
    options.SliderWidth     (1, 1) {mustBeNumeric} = 0.12;
    
    % Electrical property bounds
    options.ErBounds        (:, 2) {mustBeNumeric} = [1, 10];
    options.ErpBounds       (:, 2) {mustBeNumeric} = [-2, 1];
    options.ThkBounds       (:, 2) {mustBeNumeric} = [-1, 2];
    
    % NLayer settings
    options.GamInterp       (1, 1) {mustBeInteger, mustBePositive} = 10;
    
    % Plot settings
    options.PlotLineWidth   (1, 1) {mustBeReal} = 1.5;
    options.PlotDotWidth    (1, 1) {mustBeReal} = 1.5;
end

%% Create main figure and panels
if isfield(options, "FigureHandle")
    fig = options.FigureHandle;
else
    fig = gcf;
end
clf(fig);
fig.Position = options.MainFigDim;

plotPanel = uipanel(fig, Position=[0, 0, options.PlotPanelSize, 1]);
sliderPanel = uipanel(fig, Position=[options.PlotPanelSize, 0, 1 - options.PlotPanelSize, 1]);

erPanel = uipanel(sliderPanel, Position=[0, 0.667, 1, 0.333], Tag="er", ...
    FontSize=options.PanelFontSize, Title="Dielectric Constant (er)");

erpPanel = uipanel(sliderPanel, Position=[0, 0.333, 1, 0.333], Tag="erp", ...
    FontSize=options.PanelFontSize, Title="Dielectric Loss (erp)");

thkPanel = uipanel(sliderPanel, Position=[0, 0, 1, 0.333], Tag="thk", ...
    FontSize=options.PanelFontSize, Title="Thickness (mm)");

%% Save initial material structure
handles.f = f;
handles.NL = NL;

numLayers = size(thk, 2);

%% Create plot
ax = axes(plotPanel);

[h1, h2, h3] = zplane([], [], ax);
h1.HandleVisibility = "off";
h2.HandleVisibility = "off";
h3.HandleVisibility = "off";

xlim(ax, [-1.1, 1.1]);
ylim(ax, [-1.1, 1.1]);
hold(ax, "on");

gamPlot = cell(size(NL,2), 1);
gamFitPlot = cell(size(NL,2), 1);

% Obtain initial parameters and calculate initial values
for ii = 1:size(NL, 2)
    gamPlot{ii} = plot(ax, f{ii}, 0*f{ii}, ".", ...
        Linewidth=options.PlotLineWidth, HandleVisibility="off");
    gamFitPlot{ii} = plot(ax, f{ii}, 0*f{ii}, Linewidth=options.PlotLineWidth);
    
    gam = NL{ii}.calculate(f{ii}, er, [], thk);
    gamFit = interp(gam, options.GamInterp);
    
    gamFitPlot{ii}.XData = real(gamFit);
    gamFitPlot{ii}.YData = imag(gamFit);
    
    gamPlot{ii}.XData = real(gam);
    gamPlot{ii}.YData = imag(gam);
end

if isfield(options, "Legend")
    legend(ax, options.Legend);
end

hold(ax, "off");
handles.gamPlot = gamPlot;
handles.gamFitPlot = gamFitPlot;

%% Create sliders
createSlider = @(panel, ind, tag) uicontrol('Parent', panel, ...
    'Style', 'slider', 'Units', 'normalized', ...
    'Position', [options.SliderXPos 0.75-options.SliderYPos*(ind-1) options.SliderLength options.SliderWidth],...
    'Tag', tag);

uiSliders.erSliders = cell(numLayers, 1);
uiSliders.erpSliders = cell(numLayers, 1);
uiSliders.thkSliders = cell(numLayers, 1);
uiSliders.erRange = options.ErBounds + zeros(numLayers, 2);
uiSliders.erpRange = options.ErpBounds + zeros(numLayers, 2);
uiSliders.thkRange = options.ThkBounds + zeros(numLayers, 2);

for ii = 1:numLayers
    % Create sliders
    uiSliders.erSliders{ii} = createSlider(erPanel, ii, sprintf('er-%d',ii));
    uiSliders.erpSliders{ii} = createSlider(erpPanel, ii, sprintf('erp-%d',ii));
    uiSliders.thkSliders{ii} = createSlider(thkPanel, ii, sprintf('thk-%d',ii));
    
    % Preset initial values
    erSliderLoc = in_lerp(real(er(ii)), uiSliders.erRange(ii,:));
    if erSliderLoc > 1
        uiSliders.erSliders{ii}.Value = 1;
    elseif erSliderLoc < 0
        uiSliders.erSliders{ii}.Value = 0;
    else
        uiSliders.erSliders{ii}.Value = erSliderLoc;
    end
    
    erpSliderLoc = in_lerp(log10(abs(imag(er(ii))))/log10(10), uiSliders.erpRange(ii,:));
    if erpSliderLoc > 1
        uiSliders.erpSliders{ii}.Value = 1;
    elseif erpSliderLoc < 0
        uiSliders.erpSliders{ii}.Value = 0;
    else
        uiSliders.erpSliders{ii}.Value = erpSliderLoc;
    end
    
    thkSliderLoc = in_lerp(log10(thk(ii))/log10(10), uiSliders.thkRange(ii,:));
    if thkSliderLoc > 1
        uiSliders.thkSliders{ii}.Value = 1;
    elseif thkSliderLoc < 0
        uiSliders.thkSliders{ii}.Value = 0;
    else
        uiSliders.thkSliders{ii}.Value = thkSliderLoc;
    end
end

% Create checkbox for infinite half plane in thk panel
handles.isInfHalfPlane = uicontrol('Style', 'checkbox', 'Parent', thkPanel, ...
    'String', 'Infinite Half-Plane (Bottom Layer)', 'Units', 'normalized', ...
    'Position', [0 0.75-0.15*(numLayers) 0.5 0.1], ...
    'FontSize', 8, 'HorizontalAlignment', 'right');
handles.isInfHalfPlane.Callback = @halfPlaneValueChange;

%% Create slider labels
createTopLabel = @(panel) uicontrol('Style', 'text', 'String', 'Top Layer', ...
    'parent', panel, 'Units' , 'normalized', 'Position', [0 0.9 1 0.1], ...
    'HorizontalAlignment', 'left');


createBottomLabel = @(panel) uicontrol('Style', 'text', 'String', 'Bottom Layer', ...
    'parent', panel, 'Units' , 'normalized', 'Position', [0 0.75-options.SliderYPos*(numLayers) 0.2 0.1], ...
    'HorizontalAlignment', 'left');

createTopLabel(erPanel);
createBottomLabel(erPanel);
createTopLabel(erpPanel);
createBottomLabel(erpPanel);
createTopLabel(thkPanel);

%% Create edit field
uiEditField.erLB = cell(numLayers, 1);
uiEditField.erUB = cell(numLayers, 1);
uiEditField.erCV = cell(numLayers, 1);
uiEditField.erpLB = cell(numLayers, 1);
uiEditField.erpUB = cell(numLayers, 1);
uiEditField.erpCV = cell(numLayers, 1);
uiEditField.thkLB = cell(numLayers, 1);
uiEditField.thkUB = cell(numLayers, 1);
uiEditField.thkCV = cell(numLayers, 1);

% Create edit fields
for ii = 1:numLayers
    uiEditField.erLB{ii} = LBValueField(erPanel, ii, lerp(0, uiSliders.erRange(ii,:)));
    uiEditField.erUB{ii} = UBValueField(erPanel, ii, lerp(1, uiSliders.erRange(ii,:)));
    uiEditField.erCV{ii} = CVValueField(erPanel, ii, real(er(ii)));
    
    uiEditField.erpLB{ii} = LBValueField(erpPanel, ii, 10.^lerp(0, uiSliders.erpRange(ii,:)));
    uiEditField.erpUB{ii} = UBValueField(erpPanel, ii, 10.^lerp(1, uiSliders.erpRange(ii,:)));
    uiEditField.erpCV{ii} = CVValueField(erpPanel, ii, abs(imag(er(ii))));
    
    uiEditField.thkLB{ii} = LBValueField(thkPanel, ii, 10.^lerp(0, uiSliders.thkRange(ii,:)));
    uiEditField.thkUB{ii} = UBValueField(thkPanel, ii, 10.^lerp(1, uiSliders.thkRange(ii,:)));
    uiEditField.thkCV{ii} = CVValueField(thkPanel, ii, thk(ii));
end

handles.uiEditField = uiEditField;

%% Create slider function handle wrapper
valueChange = @(hObject, eventdata) sliderValueChanged(hObject, eventdata);

for ii = 1:numLayers
    % Set callbacks
    uiSliders.erSliders{ii}.Callback = valueChange;
    uiSliders.erpSliders{ii}.Callback = valueChange;
    uiSliders.thkSliders{ii}.Callback = valueChange;
    
    % Set listeners
    addlistener(uiSliders.erSliders{ii}, 'Value', 'PostSet', @(~, eventdata) sliderValueChanged(eventdata.AffectedObject, eventdata));
    addlistener(uiSliders.erpSliders{ii}, 'Value', 'PostSet', @(~, eventdata) sliderValueChanged(eventdata.AffectedObject, eventdata));
    addlistener(uiSliders.thkSliders{ii}, 'Value', 'PostSet', @(~, eventdata) sliderValueChanged(eventdata.AffectedObject, eventdata));
end

handles.uiSliders = uiSliders;
handles.options = options;

% Handles is an optional output parameter
if nargout == 1
    varargout{1} = handles;
end

guidata(fig, handles);
end

%% Slider value change function
function sliderValueChanged(hObject, eventdata)
handles = guidata(hObject);

% Extract value from slider and current value edit field
[er, erp, thk] = valueReader(handles);
[er_slider, erp_slider, thk_slider] = valueExtractor(handles);

panel = extractBefore(hObject.Tag,'-');
layer = str2num(extractAfter(hObject.Tag,'-'));

% Update the value being changed in the appropriate edit field
switch panel
    case 'er'
        handles.uiEditField.erCV{layer}.String = er_slider(layer);
        er(layer) = er_slider(layer);
    case 'erp'
        handles.uiEditField.erpCV{layer}.String = erp_slider(layer);
        erp(layer) = erp_slider(layer);
    case 'thk'
        handles.uiEditField.thkCV{layer}.String = thk_slider(layer);
        thk(layer) = thk_slider(layer);
end

handles = plotGam(handles, er, erp, thk);

guidata(hObject, handles);

drawnow;
end

%% Lower bound edit field creation
function LBField = LBValueField(panel, ind, initVal)
LBField = uicontrol('Style', 'edit', 'Parent', panel);
uicontrol(LBField);

pos = [0.05 0.75-0.15*(ind-1) 0.1 0.12];

LBField.Units = 'Normalized';
LBField.Position = pos;
LBField.Callback = @LBFieldChanged;
LBField.Tag = num2str(ind);
LBField.String = num2str(initVal);
end

%% Lower bound edit field callback
function LBFieldChanged(hObject, eventdata)
handles = guidata(hObject);

panel = hObject.Parent.Tag;
layer = str2num(hObject.Tag);
lowerBound = str2double(hObject.String);

if ~isnan(lowerBound)
    switch panel
        case 'er'
            currentValue = str2num(handles.uiEditField.erCV{layer}.String);
            
            if lowerBound < 1
                hObject.String = num2str(1);
                handles.uiSliders.erRange(layer,1) = 1;
            elseif lowerBound <= currentValue
                handles.uiSliders.erRange(layer,1) = lowerBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.erRange(layer,1) = currentValue;
            end
            
            guidata(hObject, handles);
            %             handles.uiSliders.erSliders{layer}.Value = in_lerp(currentValue, handles.uiSliders.erRange(layer,:));
        case 'erp'
            currentValue = str2num(handles.uiEditField.erpCV{layer}.String);
            
            if lowerBound <= currentValue
                handles.uiSliders.erpRange(layer,1) = log10(lowerBound);
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.erpRange(layer,1) = log10(currentValue);
            end
            
            guidata(hObject, handles);
            %             handles.uiSliders.erpSliders{layer}.Value = in_lerp(log10(currentValue)/log10(10), handles.uiSliders.erpRange(layer,:));
        case 'thk'
            currentValue = str2num(handles.uiEditField.thkCV{layer}.String);
            
            if lowerBound <= currentValue
                handles.uiSliders.thkRange(layer,1) = log10(lowerBound);
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.thkRange(layer,1) = log10(currentValue);
            end
            
            guidata(hObject, handles);
            %             handles.uiSliders.thkSliders{layer}.Value = in_lerp(log10(currentValue)/log10(10), handles.uiSliders.thkRange(layer,:));
    end
else
    switch panel
        case 'er'
            hObject.String = num2str(handles.uiSliders.erRange(layer,1));
        case 'erp'
            hObject.String = num2str(10.^(handles.uiSliders.erpRange(layer,1)));
        case 'thk'
            hObject.String = num2str(10.^(handles.uiSliders.thkRange(layer,1)));
    end
end

guidata(hObject, handles)
end

%% Upper bound edit field and callback
function limitField = UBValueField(panel, ind, initVal)
limitField = uicontrol('Style', 'edit', 'Parent', panel);
uicontrol(limitField);

pos = [0.745 0.75-0.15*(ind-1) 0.1 0.12];

limitField.Units = 'Normalized';
limitField.Position = pos;
limitField.Callback = @UBFieldChanged;
limitField.Tag = num2str(ind);
limitField.String = num2str(initVal);
end

function UBFieldChanged(hObject, eventdata)
handles = guidata(hObject);

panel = hObject.Parent.Tag;
layer = str2num(hObject.Tag);
upperBound = str2double(hObject.String);

if ~isnan(upperBound)
    switch panel
        case 'er'
            currentValue = str2num(handles.uiEditField.erCV{layer}.String);
            
            if upperBound >= currentValue
                handles.uiSliders.erRange(layer,2) = upperBound;
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.erRange(layer,2) = currentValue;
            end
            
            guidata(hObject, handles);
            handles.uiSliders.erSliders{layer}.Value = in_lerp(currentValue, handles.uiSliders.erRange(layer,:));
        case 'erp'
            currentValue = str2num(handles.uiEditField.erpCV{layer}.String);
            
            if upperBound >= currentValue
                handles.uiSliders.erpRange(layer,2) = log10(upperBound);
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.erpRange(layer,2) = log10(currentValue);
            end
            
            guidata(hObject, handles);
            handles.uiSliders.erpSliders{layer}.Value = in_lerp(log10(currentValue)/log10(10), handles.uiSliders.erpRange(layer,:));
        case 'thk'
            currentValue = str2num(handles.uiEditField.thkCV{layer}.String);
            
            if upperBound >= currentValue
                handles.uiSliders.thkRange(layer,2) = log10(upperBound);
            else
                hObject.String = num2str(currentValue);
                handles.uiSliders.thkRange(layer,2) = log10(currentValue);
            end
            
            guidata(hObject, handles);
            handles.uiSliders.thkSliders{layer}.Value = in_lerp(log10(currentValue)/log10(10), handles.uiSliders.thkRange(layer,:));
    end
else
    switch panel
        case 'er'
            hObject.String = num2str(handles.uiSliders.erRange(layer,2));
        case 'erp'
            hObject.String = num2str(10.^(handles.uiSliders.erpRange(layer,2)));
        case 'thk'
            hObject.String = num2str(10.^(handles.uiSliders.thkRange(layer,2)));
    end
end

guidata(hObject, handles)
end

%% Current value edit field creation
function CVField = CVValueField(panel, ind, initVal)
CVField = uicontrol('Style', 'edit', 'Parent', panel);
uicontrol(CVField);

pos = [0.87 0.75-0.15*(ind-1) 0.1 0.12];

CVField.Units = 'Normalized';
CVField.Position = pos;
CVField.Callback = @CVFieldChanged;
CVField.Tag = num2str(ind);
CVField.String = num2str(initVal);
end

%% Current value edit field callback
function CVFieldChanged(hObject, eventdata)
handles = guidata(hObject);

uiSliders = handles.uiSliders;

panel = hObject.Parent.Tag;
layer = str2num(hObject.Tag);
currentValue = str2double(hObject.String);
isOutsideBounds = 0;

if ~isnan(currentValue) || isreal(currentValue)
    [er, erp, thk] = valueReader(handles);
    numLayers = size(thk, 2);
    
    switch panel
        case 'er'
            if currentValue < 1
                hObject.String = num2str(lerp(uiSliders.erSliders{layer}.Value, uiSliders.erRange(layer,:)));
            end
            
            sliderValue = in_lerp(currentValue, uiSliders.erRange(layer,:));
            
            if sliderValue >= 0 && sliderValue <= 1
                uiSliders.erSliders{layer}.Value = sliderValue;
            elseif sliderValue < 0
                uiSliders.erSliders{layer}.Value = 0;
                isOutsideBounds = 1;
            elseif sliderValue > 1
                uiSliders.erSliders{layer}.Value = 1;
                isOutsideBounds = 1;
            else
                hObject.String = num2str(lerp(uiSliders.erSliders{layer}.Value, uiSliders.erRange(layer,:)));
            end
        case 'erp'
            sliderValue = in_lerp(log10(currentValue)/log10(10), uiSliders.erpRange(layer,:));
            
            if sliderValue >= 0 && sliderValue <= 1
                uiSliders.erpSliders{layer}.Value = sliderValue;
            elseif sliderValue < 0
                uiSliders.erpSliders{layer}.Value = 0;
                isOutsideBounds = 1;
            elseif sliderValue > 1
                uiSliders.erpSliders{layer}.Value = 1;
                isOutsideBounds = 1;
            else
                hObject.String = num2str(lerp(uiSliders.erpSliders{layer}.Value, uiSliders.erpRange(layer,:)));
            end
        case 'thk'
            sliderValue = in_lerp(log10(currentValue)/log10(10), uiSliders.thkRange(layer,:));
            
            if sliderValue >= 0 && sliderValue <= 1
                uiSliders.thkSliders{layer}.Value = sliderValue;
            elseif sliderValue < 0
                uiSliders.thkSliders{layer}.Value = 0;
                isOutsideBounds = 1;
            elseif sliderValue > 1
                uiSliders.thkSliders{layer}.Value = 1;
                isOutsideBounds = 1;
            else
                hObject.String = num2str(lerp(uiSliders.thkSliders{layer}.Value, uiSliders.thkRange(layer,:)));
            end
    end
else
    switch panel
        case 'er'
            hObject.String = num2str(lerp(uiSliders.erSliders{layer}.Value, uiSliders.erRange(layer,:)));
        case 'erp'
            hObject.String = num2str(10.^lerp(uiSliders.erpSliders{layer}.Value, uiSliders.erpRange(layer,:)));
        case 'thk'
            hObject.String = num2str(10.^lerp(uiSliders.thkSliders{layer}.Value, uiSliders.thkRange(layer,:)));
    end
end

if isOutsideBounds
    switch panel
        case 'er'
            er(layer) = currentValue;
        case 'erp'
            erp(layer) = currentValue;
        case 'thk'
            thk(layer) = currentValue;
    end
    hObject.String = currentValue;
end

handles = plotGam(handles, er, erp, thk);

handles.uiSliders = uiSliders;

guidata(hObject, handles);
end

%% Infinite half plane checkbox value change callback
function halfPlaneValueChange(hObject, eventdata)
handles = guidata(hObject);

isInfHalfPlane = handles.isInfHalfPlane.Value;

layer = size(handles.uiSliders.thkSliders, 1);

% Disable the last layer's thickness parameters textboxes and slider
if isInfHalfPlane
    handles.uiSliders.thkSliders{layer}.Enable = 'off';
    handles.uiEditField.thkCV{layer}.Enable = 'off';
    handles.uiEditField.thkLB{layer}.Enable = 'off';
    handles.uiEditField.thkUB{layer}.Enable = 'off';
else
    handles.uiSliders.thkSliders{layer}.Enable = 'on';
    handles.uiEditField.thkCV{layer}.Enable = 'on';
    handles.uiEditField.thkLB{layer}.Enable = 'on';
    handles.uiEditField.thkUB{layer}.Enable = 'on';
end

[er, erp, thk] = valueExtractor(handles);

panel = extractBefore(hObject.Tag,'-');
layer = str2num(extractAfter(hObject.Tag,'-'));

switch panel
    case 'er'
        handles.uiEditField.erCV{layer}.String = er(layer);
    case 'erp'
        handles.uiEditField.erpCV{layer}.String = erp(layer);
    case 'thk'
        handles.uiEditField.thkCV{layer}.String = thk(layer);
end

handles = plotGam(handles, er, erp, thk);


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
v_arr = range(:,1) + m.*(range(:,end) - range(:,1));
end

%% Value extractor
function [er, erp, thk] = valueExtractor(handles)
valueExtracted = @(s) s.Value;
er = lerp_arr(cellfun(valueExtracted, handles.uiSliders.erSliders), handles.uiSliders.erRange);
erp = 10.^lerp_arr(cellfun(valueExtracted, handles.uiSliders.erpSliders), handles.uiSliders.erpRange);
thk = 10.^lerp_arr(cellfun(valueExtracted, handles.uiSliders.thkSliders), handles.uiSliders.thkRange);
end

%% Value reader from current value edit field
function [er, erp, thk] = valueReader(handles)
% Read current value from edit field which is string stored in cell array
valueRead = @(s) s.String;
er_str = cellfun(valueRead, handles.uiEditField.erCV, UniformOutput=false);
erp_str = cellfun(valueRead, handles.uiEditField.erpCV, UniformOutput=false);
thk_str = cellfun(valueRead, handles.uiEditField.thkCV, UniformOutput=false);

% Convert string to number
er = cellfun(@str2num, er_str);
erp = cellfun(@str2num, erp_str);
thk = cellfun(@str2num, thk_str);
end

%% Calculate and plot S-Parameters
function handles = plotGam(handles, er, erp, thk)
er_in = er - 1j*erp;
thk_in = thk;

if handles.isInfHalfPlane.Value
    thk_in(end) = inf;
end

for ii = 1:size(handles.NL, 2)
    gam = handles.NL{ii}.calculate(handles.f{ii}, er_in.', [], thk_in.');
    gamFit = interp(gam, 10);
    
    handles.gamFitPlot{ii}.XData = real(gamFit);
    handles.gamFitPlot{ii}.YData = imag(gamFit);
    
    handles.gamPlot{ii}.XData = real(gam);
    handles.gamPlot{ii}.YData = imag(gam);
end
end