clc;
clear;
close all;

%% Inputs
a = 2;
b = 1;

x(:, 1) = a * linspace(0, 1, 2001);
y(1, :) = b * linspace(0, 1, 1001);

m(1, 1, :) = [1, 2, 3, 4, 1, 2, 3, 4];
n(1, 1, :) = [0, 0, 0, 0, 1, 1, 1, 1];

dim4 = @(x) reshape(x, size(x, 1), size(x, 2), 1, []);

%% Calculate Fields Everywhere
HxTE = (m ./ a) .* sin(m .* pi .* x ./ a) .* cos(n .* pi .* y ./ b);
HyTE = (n ./ b) .* cos(m .* pi .* x ./ a) .* sin(n .* pi .* y ./ b);
ExTE = HyTE;
EyTE = HxTE;

ExTM = (m ./ a) .* cos(m .* pi .* x ./ a) .* sin(n .* pi .* y ./ b);
EyTM = (n ./ b) .* sin(m .* pi .* x ./ a) .* cos(n .* pi .* y ./ b);
HxTM = EyTM;
HyTM = ExTM;

%% Integrate
mat1 = innerProduct(HxTE, dim4(HxTE), [1, 2], SummationMode="Mean") ...
    +  innerProduct(HyTE, dim4(HyTE), [1, 2], SummationMode="Mean");
abs(squeeze(mat1))

mat2 = innerProduct(HyTE, dim4(HxTE), [1, 2], SummationMode="Mean") ...
    -  innerProduct(HxTE, dim4(HyTE), [1, 2], SummationMode="Mean");
abs(squeeze(mat2))


%% Plot
% figure;
% showImage(x, y, ExTE(:, :, 5), DisplayFormat="Magnitude");

