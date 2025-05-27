function [] = setupSliders(fig, sliderPanels, er, ur, thk, ...
    erpRange, erTanRange, urpRange, urRanRange, thkRange, ...
    structureUpdateFun)
%Setup sliders for nLayerViewer.
% This is for nLayerViewer.
% 
% Author: Matt Dvorsky

arguments
    fig(1, 1) matlab.ui.Figure;
    sliderPanels(:, 1);
    er;
    ur;
    thk;
    erpRange   (:, 2);
    erTanRange (:, 2);
    urpRange   (:, 2);
    urRanRange (:, 2);
    thkRange   (:, 2);
    structureUpdateFun;
end

%% Get Current, Min, and Max Values
initialStructureArray = structureToArray(er, ur, thk);

zeroArr = zeros(size(initialStructureArray, 1), 1);
minValArray = [...
    zeroArr + erpRange(:, 1), ...
    zeroArr + erTanRange(:, 1), ...
    zeroArr + urpRange(:, 1), ...
    zeroArr + urRanRange(:, 1), ...
    zeroArr + thkRange(:, 1), ...
    ];
maxValArray = [...
    zeroArr + erpRange(:, 2), ...
    zeroArr + erTanRange(:, 2), ...
    zeroArr + urpRange(:, 2), ...
    zeroArr + urRanRange(:, 2), ...
    zeroArr + thkRange(:, 2), ...
    ];

%% Create Sliders
for ss = 1:5
    for ii = 1:numel(thk)
        addSlider(sliderPanels{ss}, ...
            initialStructureArray(ii, ss), ...
            [minValArray(ii, ss), maxValArray(ii, ss)], ...
            SliderHeightOffset=(ii - 0.5)*0.2, ...
            SliderHeight=0.15, ...
            ValueChangedHandler={@sliderHandler, ...
                fig, structureUpdateFun, ii, ss});
    end
end

end




%% Helper Functions
function sliderHandler(val, fig, structureUpdateFun, layerInd, sInd)
    fig.UserData.structure(layerInd, sInd) = val;

    try
        [er, ur, thk] = arrayToStructure(fig.UserData.structure);
        structureUpdateFun(fig, er, ur, thk);
    catch ex
        panel = fig.Children(end);
        oldColor = panel.BackgroundColor;
        for ii = 1:2
            panel.BackgroundColor = "red";
            pause(0.1);
            panel.BackgroundColor = oldColor;
            pause(0.1);
        end

        rethrow(ex);
    end
end
