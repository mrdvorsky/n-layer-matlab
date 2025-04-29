clc;
clear;
close all;

%% Inputs
kc = 1;
kr(1, :) = linspace(0.01, 20, 100);

%% Nodes
[x, w] = fejer2_halfOpen(100000, 1);
% [x, w] = fejer2(1000000, 0, 1000);


fun = @(x, kr) besselj(0, kr.*x) .* x ./ (x.^2 - kc.^2);

val = sum(w .* fun(x, kr), 1);


for ii = 1:numel(kr)
    val(ii) = integral(@(x) fun(x, kr(ii)), 0, 1j);
end

figImg = figure;
showImage(1:10, 1:10, rand(10));
while true
    [x, y, v] = impixel();

    figConv = figure;
    showImage(1:10, 1:10, rand(10));

    button = waitForSpecificButtonPress('c', 'q');
    if button == 'c'
        close(figConv);
        continue;
    elseif button == 'q'
        break;
    end
end


function [button] = waitForSpecificButtonPress(varargin)
    fig = gcf;
    while true
        if waitforbuttonpress() == 1    % If keyboard
            if any(strcmp(string(varargin), fig.CurrentCharacter))
                button = fig.CurrentCharacter;
                return;
            end
        end
    end
end


% plot(kr, val);
% ylim([-1, 5])


