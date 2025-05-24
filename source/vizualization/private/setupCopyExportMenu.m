function [] = setupCopyExportMenu(fig, plotAxis, structureAxis)
%Sets up the menu bar for copying and exporting the figure.
% This is for nLayerViewer.
% 
% Author: Matt Dvorsky

arguments
    fig(1, 1) matlab.ui.Figure;
    plotAxis(1, 1) matlab.graphics.axis.Axes;
    structureAxis(1, 1) matlab.graphics.axis.Axes;
end

%% Create Menu Bar Entry
copyMenu = uimenu(fig, Text="Copy/Export Figure");

%% Add Items
uimenu(copyMenu, ...
    Text="Copy Polar Plot", ...
    MenuSelectedFcn={@copyPolarPlot, plotAxis});
uimenu(copyMenu, ...
    Text="Copy Structure Definition", ...
    MenuSelectedFcn={@copyStructure, structureAxis});
uimenu(copyMenu, ...
    Text="Copy Both", ...
    MenuSelectedFcn={@copyPolarPlotAndStructure, plotAxis});

uimenu(copyMenu, ...
    Text="Export Polar Plot", ...
    Separator="on", ...
    MenuSelectedFcn={@exportPolarPlot, plotAxis});
uimenu(copyMenu, ...
    Text="Export Structure Definition", ...
    MenuSelectedFcn={@exportStructure, structureAxis});
uimenu(copyMenu, ...
    Text="Export Both", ...
    MenuSelectedFcn={@exportPolarPlotAndStructure, plotAxis, structureAxis});

end



%% Handler Functions
function copyPolarPlot(~, ~, plotAxis)
    copygraphics(plotAxis, ...
        Resolution=200, ...
        ContentType="image");
end

function copyStructure(~, ~, structureAxis)
    copygraphics(structureAxis, ...
        Resolution=200, ...
        ContentType="image");
end

function copyPolarPlotAndStructure(~, ~, plotAxis)
    copygraphics(plotAxis.Parent, ...
        Resolution=200, ...
        ContentType="image");
end


function exportPolarPlot(~, ~, plotAxis)
    newFig = figure();
    newPlotAxis = copyobj(plotAxis, newFig);
    newPlotAxis.Position = [0.1300 0.1100 0.7750 0.8150];

    if isprop(plotAxis, "Legend")
        legend(newPlotAxis);
    end
end

function exportStructure(~, ~, structureAxis)
    newFig = figure();
    newStructureAxis = copyobj(structureAxis, newFig);
    newStructureAxis.Children(1).FontSize = 9;
    newStructureAxis.Position = [0.1, 0.1, 0.8, 0.8];
end

function exportPolarPlotAndStructure(~, ~, plotAxis, structureAxis)
    plotPanel = plotAxis.Parent;
    plotPanel.Units = "pixels";

    newFig = figure(Position=plotPanel.Position);

    newPlotAxis = copyobj(plotAxis, newFig);
    newStructureAxis = copyobj(structureAxis, newFig);
    newStructureAxis.Children(1).FontSize = 9;

    if isprop(plotAxis, "Legend")
        legend(newPlotAxis);
    end
end

