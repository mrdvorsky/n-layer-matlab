clc;
clear;
close all;

%% Inputs
er = {4 - 0.4j};
ur = {1};
thk = {1.2};

f = linspace(0.1, 18, 201);

numModes = 5;

%% nLayerCircular
tic;
NL1 = nLayerCoaxial(0, numModes - 1, waveguideBand="N");
% NL1.waveguideRo = 7;
% NL1.waveguideRi = 2;
NL1.calculate(f(1), er, ur, thk);
toc;

%% Mode Plotting
% for ii = 1:numel(NL1.modeStructs)
%     nLayer.plotModeStruct(NL1.modeStructs(ii));
% end

%% Calculate
tic;
gam1 = NL1.calculate(f, er, ur, thk);
toc;

%% Plot
figure;
plot(gam1, "-", LineWidth=1.5);
hold on;
zplane([]);
grid on;
% legend(["Original", "New"]);


% nLayer.plotModeStruct(NL2.modeStructs);

% figure;
% nLayerViewer(4, 1, 2, NL2, f);

