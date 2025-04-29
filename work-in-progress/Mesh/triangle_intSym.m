clc;
clear;
close all;

%% Inputs
syms x y kx ky;

%% Integration
int((1 - y) .* exp(-1j .* kx .* x) .* exp(-1j * ky .* y), y, 1)




