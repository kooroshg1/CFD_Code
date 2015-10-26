clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
CSA_dTdX = dlmread('CSA_dTdX.txt');
csa_dTdX = CSA_dTdX;
csa_dTdX(12:20, 12:20) = csa_dTdX(12:20, 12:20) + 2 * (rand(9,9) - 0.5);

dlmwrite('CSA_dTdX_Old.txt', csa_dTdX)

X = dlmread('complex_step_X.txt');
Y = dlmread('complex_step_Y.txt');

figure,
plot(X(:, 16), CSA_dTdX(:, 16), 'ro',...
     X(:, 16), csa_dTdX(:, 16), 'k')
legend('CSA', 'CS')
