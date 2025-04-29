clc;
clear;
close all;

%% Inputs
orderNum = 4;

%% Get Shape Functions
shapeFun = fem_createShapeFunctions(orderNum);

%% Get FEM Matrices
[T, Q] = fem_generateElementMatrices(orderNum);

