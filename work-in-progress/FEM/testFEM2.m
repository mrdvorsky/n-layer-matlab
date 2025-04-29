clc;
clear;
close all;

%% Inputs
model = createpde();

d1 = decsg([1, 0, 0, 1].');
geometryFromEdges(model, d1);
applyBoundaryCondition(model, "dirichlet", Edge=1:4, u=0);
specifyCoefficients(model, c=1, d=1, a=0, f=0, m=0);

generateMesh(model, "GeometricOrder", "quadratic");


results = solvepdeeig(model, [0, 50]);

%% Plotting
for ii = 1:numel(results.Eigenvalues)
    figure;
    pdeplot(model, XYData=results.Eigenvectors(:, ii), ZData=results.Eigenvectors(:, ii));
end

