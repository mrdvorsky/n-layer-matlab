function [modes] = setModes(O, maxM, maxN)
%SETMODES Update modes in nLayerRectangular object and return a list

O.modes = [reshape((1:2:maxM).' + 0*(0:2:maxN), [], 1), ...
    reshape(0*(1:2:maxM).' + (0:2:maxN), [], 1)];
modes = O.modes;

end


