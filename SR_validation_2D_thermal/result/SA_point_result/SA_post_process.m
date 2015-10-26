clc;
clear all;
close all
% ----------------------------------------------------------------------- %
% T0 = dlmread('T0.txt');
% T1 = dlmread('T1.txt');
% X = dlmread('X.txt');
% Y = dlmread('Y.txt');
% 
% % dTdX = (T1 - T0) / 0.001;
% dTdX = imag(T1) / 0.001;
% dlmwrite('complex_step_dTdX.txt', dTdX);
% dlmwrite('complex_step_X.txt', X);
% dlmwrite('complex_step_Y.txt', Y);

CS_dTdX = dlmread('complex_step_dTdX.txt');
SA_dTdX = dlmread('SA_dTdX.txt');
X = dlmread('complex_step_X.txt');
Y = dlmread('complex_step_Y.txt');

figure,
plot(X(:, 26), SA_dTdX(:, 26), 'ro',...
     X(:, 26), CS_dTdX(:, 26), 'k')
