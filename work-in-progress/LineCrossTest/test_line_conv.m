clc;
clear;
close all;

%% Inputs
p1_1 = [10, 10];
p1_2 = [20, 20];

p2_1 = [5, 1];
p2_2 = [10, -1];

line1_x = linspace(p1_1(1), p1_2(1), 100);
line1_y = linspace(p1_1(2), p1_2(2), 100);

line2_x = linspace(p2_1(1), p2_2(1), 100);
line2_y = linspace(p2_1(2), p2_2(2), 100);

%% Plot
figure;
plot(line1_x, line1_y, LineWidth=1.5);
hold on;
plot(line2_x, line2_y, LineWidth=1.5);
grid on;





