function [] = nLayerViewer(er, ur, thk, NL, f, options)
%Creates an interactive plot for visualizing multi-layered structure measurements.
% Creates a figure with an S-parameter plot, as well as sliders to change
% the multi-layered structure paramaters.
%
% ===== Basic Usage =====
%   figure;
%   nLayerViewer(er, ur, thk, NL, f);
%
% ===== Multiple nLayerForward Objects =====
%   figure;
%   nLayerViewer(er, ur, thk, NL1, f1, NL2, f2, ...);
%
%
% Inputs:
%   er - Vector of default complex permittivity for each layer.
%   ur - Same as "er", but for complex permeability.
%   thk - Same as "er", but for layer thicknesses.
%   NL (Repeating) - The nLayerForward solver object(s) to use.
%   f (Repeating) - Vector(s) of frequencies. Only the min and max values
%       will be used to calculate the plot frequencies.
%
% Outputs:
%   (none)
%
% Options (name-value pairs):
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
% Author: Matt Dvorsky

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
    options.Figure     (1, 1) matlab.ui.Figure;
    options.FigureSize (1, 4) {mustBeReal} = [100, 100, 1000, 600];

    % Legend
    options.ShowLegend      (1, 1) logical = true;
    options.PlotPanelWidth  (1, 1) {mustBeFractional} = 0.6;
    options.PanelFontSize   (1, 1) {mustBePositive} = 10;

    % Slider parameters
    options.SliderXPos      (1, 1) {mustBeFractional} = 0.15;
    options.SliderYPos      (1, 1) {mustBeFractional} = 0.15;
    options.SliderLength    (1, 1) {mustBeFractional} = 0.6;
    options.SliderWidth     (1, 1) {mustBeFractional} = 0.12;

    % Structure property bounds
    options.ErpBounds       (:, 2) {mustBeReal, mustBeFinite} = [1, 10];
    options.ErLossTanBounds (:, 2) {mustBeReal, mustBeFinite} = [0, 1];
    options.UrpBounds       (:, 2) {mustBeReal, mustBeFinite} = [1, 10];
    options.UrLossTanBounds (:, 2) {mustBeReal, mustBeFinite} = [0, 1];
    options.ThkBounds (:, 2) {mustBeNonnegative, mustBeFinite} = [0, 10];

    % Frequency Settings
    options.NumFrequencies      (1, 1) {mustBeInteger, mustBePositive} = 801;
    options.NumFrequencyMarkers (1, 1) {mustBeInteger, mustBePositive} = 21;

    % Plot settings
    options.PlotLineWidth   (1, 1) {mustBePositive} = 1;
    options.PlotMarkerSize  (1, 1) {mustBePositive} = 5;
    options.PlotMarkerType  (1, 1) string = ".";
end
mustHaveAtLeastOneRepeatingArg(NL);

%% Check Inputs
[er, ur, thk] = nLayer.validateStructure(er, ur, thk, ...
    RequireScalarValuesPerLayer=false, ...
    CheckStructureValues=false);

if isfield(options, "Figure")
    fig = options.Figure;
else
    fig = gcf();
end

%% Setup Figure, Panels, etc.
[plotPanel, sliderPanels] = setupFigureAndPanels(fig, ...
    options.FigureSize, ...
    options.PlotPanelWidth, ...
    options.PanelFontSize);
setupCopyExportMenu(fig);
[plotAxis, structureTextHandle] = setupPlotPanel(plotPanel, numel(thk));
setupSliders(fig, sliderPanels, er, ur, thk, ...
    options.ErpBounds, ...
    options.ErLossTanBounds, ...
    options.UrpBounds, ...
    options.UrLossTanBounds, ...
    options.ThkBounds, ...
    @structureUpdateFun);

%% Setup Initial Polar Plot
numFreqsPerMarker = ceil(options.NumFrequencies ./ options.NumFrequencyMarkers);
numF =  numFreqsPerMarker .* options.NumFrequencyMarkers;
fMarkerInds = (1:options.NumFrequencyMarkers) .* numFreqsPerMarker;

plotUpdateFuns = cell(numel(NL), 1);
hold(plotAxis, "on");
for ii = 1:numel(NL)
    plotLabels = NL{ii}.getOutputLabels();
    f{ii} = linspace(min(f{ii}), max(f{ii}), numF);

    [lineHandles, plotUpdateFuns{ii}] = plotComplex(...
        f{ii}, ...
        zeros(numel(f{ii}, numel(plotLabels))), ...
        "", ...
        AddCustomDataTips=true, ...
        Axis=plotAxis, ...
        Linewidth=options.PlotLineWidth, ...
        Marker=options.PlotMarkerType, ...
        MarkerSize=options.PlotMarkerSize, ...
        MarkerIndices=fMarkerInds);

    for pp = 1:numel(lineHandles)
        lineHandles(pp).DisplayName = plotLabels(pp);
    end
end

if options.ShowLegend
    legend(plotAxis);
end

datacursormode(fig, "on");

%% Store Data in Figure
fig.UserData.structure = structureToArray(er, ur, thk);
fig.UserData.NL = cellfun(@copy, NL, UniformOutput=false);
fig.UserData.f = f;
fig.UserData.plotUpdateFuns = plotUpdateFuns;
fig.UserData.structureTextHandle = structureTextHandle;

%% Run Updater Once
structureUpdateFun(fig, er, ur, thk);

end




% function figureResizeCallback(fig, ~)
%     handles = guidata(fig);
%     figMinRatio = min(fig.Position([3, 4]) ./ [1000, 600]);
%     handles.structureText.FontSize = round(9 .* figMinRatio);
% end

%% Update Functions
function structureUpdateFun(fig, er, ur, thk)
    NL = fig.UserData.NL;
    f = fig.UserData.f;
    plotUpdateFuns = fig.UserData.plotUpdateFuns;

    for ii = 1:numel(NL)
        gam = NL{ii}.calculate(f{ii}, er, ur, thk);
        plotUpdateFuns{ii}(gam(:, :));
    end

    [~, fig.UserData.structureTextHandle.String] = ...
        nLayer.printStructure(er, ur, thk, Title="");
end
