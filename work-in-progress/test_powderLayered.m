clc;
clear;
close all;

%% Inputs
f(:, 1) = linspace(26.5, 40, 201);

numLayers = 61;

er1 = 15 - 2.5j;
er2 = 1;

ur1 = 1;
ur2 = 0.45 - 0.05j;

thk1 = 0.1;
thk2 = 0.2;

erGuess = 6 - 0.6j;
urGuess = 0.7 - 0.07j;

%% Make Structure
er = repmat({er1}, numLayers, 1);
er(2:2:end) = repmat({er2}, floor(0.5*numLayers), 1);

ur = repmat({ur1}, numLayers, 1);
ur(2:2:end) = repmat({ur2}, floor(0.5*numLayers), 1);

thk = repmat({thk1}, numLayers, 1);
thk(2:2:end) = repmat({thk2}, floor(0.5*numLayers), 1);

nLayer.printStructure(er, ur, thk);

%% nLayer
NL = nLayerFilledRectangular(1, 0, waveguideBand="Ka");
NL.checkStructureValues = false;
gamMeas = NL.calculate(f, er, ur, thk);

%% Inverse
NLsolver = nLayerInverse(1, verbosity=1);
NLsolver.setLayersToSolve(Erp=[1], Erpp=[1], Urp=[1], Urpp=[1]);
NLsolver.rangeMin_urp = [0.1];
NLsolver.setInitialValues(Er=erGuess, Ur=urGuess, Thk=sum(cell2mat(thk)));

NLsolver.printStructureParameters();

[Params, ~, Uncert] = NLsolver.solveStructure(NL, f, gamMeas);

NLsolver.printStructureParameters(Params, Uncert);

%% Plotting
figure;
nLayerViewer(Params.er, Params.ur, Params.thk, NL, f);
hold on;
h = plot(gamMeas(:, :), ":", Linewidth=1.5);
set(h, {'DisplayName'}, cellstr(compose("%s (Meas)", NL.getOutputLabels().')));


