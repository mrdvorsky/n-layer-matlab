clc;
clear;
close all;

%% Inputs
f = [6.22];

%% Define Structure
soil_er = 3.35 * (1 - 0.02j);
er = [1, (soil_er), 11.2 - 1.5j, soil_er];
thk = [0, 39, 9, 10000];

%% Create nLayer Object
NL = nLayerPlaneWave(verbosity=1, outputIndices=1);

%% nLayer Inverse
NL.printStructure(er, [], [thk(1:end-1), inf]);
gam0 = NL.calculate(f, er, [], thk);

NLsolver = nLayerInverse(length(thk), Verbosity=0);
NLsolver.setInitialValues(ErValue=real(er), ErpValue=-imag(er), ThkValue=thk);
NLsolver.setLayersToSolve(ErLayers=[], ErpLayers=[], ThkLayers=[2, 3]);

erpGuesses = [11.2, linspace(10, 15, 11)];
for ii = 1:numel(erpGuesses)
    NLsolver.erInitialValue(3) = erpGuesses(ii);
    [erOut, ~, thkOut, gamOut] = NLsolver.solveStructure(NL, f, gam0);
%     abs(gam0 - gamOut)

    fprintf("(er3' = %.1f: thk2 = %.5f, thk3 = %.5f) -> gam = %.4f - j%.4f\n", ...
        real(erOut(3)), thkOut(2), thkOut(3), real(gamOut), -imag(gamOut));
end


return;

%% Plot Plane 1
erpp_1(:, 1) = linspace(0, 3, 101);
erp_1(1, :) = linspace(1, 15, 201);

gamError1 = 0*erpp_1 + 0*erp_1;
for ii = 1:numel(erpp_1)
    for jj = 1:numel(erp_1)
        er1 = er;
        thk1 = thk;
        er1(2) = complex(erp_1(jj), -erpp_1(ii));
        gamError1(ii, jj) = NL.calculate(f, er1, [], thk1) - gam0;
    end
end

figure;
showImage(erpp_1, erp_1, gamError1, DisplayFormat="dB");
xlabel("er''");
ylabel("er'");
axis normal;

%% Plot Plane 2
d3_2(:, 1) = linspace(0, 20, 101);
d2_2(1, :) = linspace(0, 40, 201);

gamError2 = 0*d3_2 + 0*d2_2;
for ii = 1:numel(d3_2)
    for jj = 1:numel(d2_2)
        er1 = er;
        thk1 = thk;
        thk1(1) = d2_2(jj);
        thk1(2) = d2_2(ii);
        gamError2(ii, jj) = NL.calculate(f, er1, [], thk1) - gam0;
    end
end

figure;
showImage(d3_2, d2_2, gamError2, DisplayFormat="dB");
xlabel("d3");
ylabel("d2");
axis normal;

%% Plot Plane 3
erp_3(:, 1) = linspace(1, 15, 101);
d3_3(1, :) = linspace(0, 20, 201);

gamError3 = 0*erpp_1 + 0*erp_1;
for ii = 1:numel(erp_3)
    for jj = 1:numel(d3_3)
        er1 = er;
        thk1 = thk;
        er1(2) = complex(erp_3(ii), imag(er1(2)));
        thk1(2) = d3_3(jj);
        gamError3(ii, jj) = NL.calculate(f, er1, [], thk1) - gam0;
    end
end

figure;
showImage(erp_3, d3_3, gamError3, DisplayFormat="dB");
xlabel("d3");
ylabel("d2");
axis normal;
