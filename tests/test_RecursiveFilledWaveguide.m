clc
clear
close all

%% Script to test recursive filled rectangular to old filled rectangular
% Randomly generates structures with varying number of layers layers, 
% permittivities, and permeatilities. Then these structures are evaluated 
% with both models and the maximum error is printed out at the end. 

% *** NOTE *** 
% There appears to be a bug in the old formulatoins which does not
% correctly handle magnetic layers. When fixed in the future uncomment "ur" 
% in the for loop and verify good match. 

NLFnew = nLayerFilledRectangular(1,0,waveguideBand="Ka");
NLFold = nLayerFilledRectangular_old(1,0);
NLFold.waveguideA = NLFnew.waveguideA;
NLFold.waveguideB = NLFnew.waveguideB;

f = linspace(32,40,1001);
tand = 0.1;
numLayers = 5;
numRuns = 100;
ur = ones([1,numLayers]);
error = zeros([numRuns,1]);

for ii=1:(numRuns)
    er = randi([1,20],[numLayers,1]) - 1j*((1e-4-1).*rand(numLayers,1) + 1);
    % ur = randi([1,20],[numLayers,1]) - 1j*((1e-4-1).*rand(numLayers,1) + 1);
    thk = (0.1-30).*rand(numLayers,1) + 30;
    gamOld = NLFold.calculate(f,er,ur,thk);
    gamNew = NLFnew.calculate(f,er,ur,thk);
    error(ii) = rmse(gamOld,gamNew,'all');
    fprintf("Error=%g\n", error(ii));
end

fprintf("Maximum Error = %g dB",db(max(error,[],'all')));

