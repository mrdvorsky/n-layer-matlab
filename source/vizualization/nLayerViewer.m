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
%   Figure (gcf) - Handle to existing figure to use.
%   FigureSize ([0.15, 0.15, 0.7, 0.7]) - Figure "Position" field to set,
%       in normalized units.
%   FontSize (10) - Font size to use.
%
%   ErpBounds       ([1, 10]) - Permittivity real part bounds.
%   ErLossTanBounds ([0, 1])  - Permittivity loss tangent bounds.
%   UrpBounds       ([1, 10]) - Permeability real part bounds.
%   UrLossTanBounds ([0, 1])  - Permeability loss tangent bounds.
%   ThkBounds       ([0, 10)  - Thickness bounds.
%
%   PlotPanelWidth (0.6) - Relative width of plot panel.
%   PlotLineWidth  (1)   - Line width for polar plot.
%   PlotMarkerSize (8)   - Marker size for polar plot.
%   PlotMarkerType (".") - Type of marker for polar plot.
%   ShowLegend    (true) - Whether to show plot legend.
%
%   NumFrequencies      (801) - Total number of frequency points to use.
%   NumFrequencyMarkers (31)  - Number of marker points on plot.
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
    options.Figure     (1, 1) matlab.ui.Figure;
    options.FigureSize (1, 4) {mustBeFractional} = [0.15, 0.15, 0.7, 0.7];
    options.FontSize   (1, 1) {mustBePositive} = 10;

    options.ErpBounds       (:, 2) {mustBeReal, mustBeFinite} = [1, 10];
    options.ErLossTanBounds (:, 2) {mustBeReal, mustBeFinite} = [0, 1];
    options.UrpBounds       (:, 2) {mustBeReal, mustBeFinite} = [1, 10];
    options.UrLossTanBounds (:, 2) {mustBeReal, mustBeFinite} = [0, 1];
    options.ThkBounds (:, 2) {mustBeNonnegative, mustBeFinite} = [0, 10];

    options.NumFrequencies      (1, 1) {mustBeInteger, mustBePositive} = 801;
    options.NumFrequencyMarkers (1, 1) {mustBeInteger, mustBePositive} = 31;

    options.PlotPanelWidth  (1, 1) {mustBeFractional} = 0.6;
    options.PlotLineWidth   (1, 1) {mustBePositive} = 1;
    options.PlotMarkerSize  (1, 1) {mustBePositive} = 8;
    options.PlotMarkerType  (1, 1) string = ".";
    options.ShowLegend      (1, 1) logical = true;
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
    options.FontSize);

[plotAxis, structureTextHandle] = setupPlotPanel(plotPanel, numel(thk));

setupCopyExportMenu(fig, plotAxis, structureTextHandle.Parent);

setupSliders(fig, sliderPanels, er, ur, thk, ...
    options.ErpBounds, ...
    options.ErLossTanBounds, ...
    options.UrpBounds, ...
    options.UrLossTanBounds, ...
    options.ThkBounds, ...
    @structureUpdateFun);

%% Setup Initial Polar Plot
if options.NumFrequencyMarkers >= 2
    numFreqsPerMarker = ceil((options.NumFrequencies - 1) ...
        ./ (options.NumFrequencyMarkers - 1));
    numF = 1 + numFreqsPerMarker.*(options.NumFrequencyMarkers - 1);
    fMarkerInds = 1 + (0:options.NumFrequencyMarkers - 1).*(numFreqsPerMarker);
else
    numF = options.NumFrequencies;
    fMarkerInds = [];
end

plotUpdateFuns = cell(numel(NL), 1);
hold(plotAxis, "on");
for ii = 1:numel(NL)
    plotLabels = NL{ii}.getOutputLabels();
    f{ii} = linspace(min(f{ii}), max(f{ii}), numF);

    [lineHandles, plotUpdateFuns{ii}] = plotComplex(...
        f{ii}, ...
        zeros(numel(f{ii}), numel(plotLabels)), ...
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
