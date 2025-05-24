function [] = setupCopyExportMenu(fig)
%Sets up the menu bar for copying and exporting the figure.
% This is for nLayerViewer.
% 
% Author: Matt Dvorsky

arguments
    fig(1, 1) matlab.ui.Figure;
end

%% Create Menu Bar Entry
copyMenu = uimenu(fig, Text="Copy/Export Figure");

%% Add Items
uimenu(copyMenu, ...
    Text="Copy Polar Plot", ...
    MenuSelectedFcn=@copyFigure);
uimenu(copyMenu, ...
    Text="Copy Structure Definition", ...
    MenuSelectedFcn=@copyStructure);
uimenu(copyMenu, ...
    Text="Copy Both", ...
    MenuSelectedFcn=@copyFigureAndStructure);

uimenu(copyMenu, ...
    Text="Export Polar Plot", ...
    MenuSelectedFcn=@exportFigure, Separator="on");
uimenu(copyMenu, ...
    Text="Export Structure Definition", ...
    MenuSelectedFcn=@exportStructure);
uimenu(copyMenu, ...
    Text="Export Both", ...
    MenuSelectedFcn=@exportFigureAndStructure);

end



%% Handler Functions
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

