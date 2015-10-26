clc;
% clear all;
close all;
% ------------------------------------------------------------ %
X = dlmread('X.txt');
Y = dlmread('Y.txt');
T = dlmread('T.txt');

loc = size(X,1)/2;
x1 = sqrt(0.1^2 - (Y(2,loc) -0.5).^2) + 0.5
x2 = -sqrt(0.1^2 - (Y(2,loc) -0.5).^2) + 0.5

figure,
plot(X(:,loc), T(:,loc), 'bo',...
       [x1, x2], [1, 1], 'r+')
   
figure,
contourf(X,Y,T)

% xloc = 0.5;
% yloc = 0.5;
% 
% figure,
% plot(X(:,size(X,1)/2), T(:,size(T,1)/2), 'bo',...
%      xloc, 1, 'r+')
%    
% figure,
% contourf(X,Y,T)
