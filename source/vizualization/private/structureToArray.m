function [structureArray] = structureToArray(er, ur, thk)
%Convert structure from "structure form" to 2D array form.
% This is for nLayerViewer.
% 
% Author: Matt Dvorsky

arguments
    er;
    ur;
    thk;
end

%% Convert
erpArr  = real(cell2mat(er(:)));
erppArr = imag(cell2mat(er(:)));
urpArr  = real(cell2mat(ur(:)));
urppArr = imag(cell2mat(ur(:)));
thkArr  = cell2mat(thk(:));

structureArray = [...
    erpArr, ...
    -erppArr ./ erpArr, ...
    urpArr, ...
    -urppArr ./ urpArr, ...
    thkArr, ...
    ];

end
