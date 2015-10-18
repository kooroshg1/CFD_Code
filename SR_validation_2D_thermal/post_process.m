clc;
% clear all;
close all;
% ------------------------------------------------------------ %
X = dlmread('X.txt');
Y = dlmread('Y.txt');
T = dlmread('T.txt');

loc = 19;
x1 = sqrt(0.1^2 - (Y(2,loc) -0.5).^2) + 0.5
x2 = -sqrt(0.1^2 - (Y(2,loc) -0.5).^2) + 0.5

figure,
plot(X(:,loc), T(:,loc),...
       [x1, x2], [1, 1], 'ro')
