clc;
clear;
close all;

%% Inputs
NL = nLayerCircular(3, 3, waveguideBand="Ka_TE01", ...
    modeSymmetryX="None", ...
    modeSymmetryY="None", ...
    modeSymmetryAxial="None");

for ii = 1:NL.numModes
    NL.waveguideModes(ii).validateModeSymmetry();

    figure;
    NL.waveguideModes(ii).showMode();
    title(sprintf("x='%s', y='%s', axial='%s'", ...
        NL.waveguideModes(ii).symmetryX, ...
        NL.waveguideModes(ii).symmetryY, ...
        NL.waveguideModes(ii).symmetryAxial));
end


