clc;
clear all;
close all
% ----------------------------------------------------------------------- %
fontsize = 40;
linewdith = 5;
markersize = 15;
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
% 
% figure,
% surf(X, Y, dTdX)

CS_dTdX = dlmread('complex_step_dTdX.txt');
CSA_dTdX = dlmread('CSA_dTdX.txt');
CSA_dTdX_Old = dlmread('CSA_dTdX_Old.txt');
X = dlmread('complex_step_X.txt');
Y = dlmread('complex_step_Y.txt');

figure,
plot(X(:, 16), CSA_dTdX(:, 16), 'ro',...
     X(:, 16), CSA_dTdX_Old(:, 16), 'b+',...
     X(:, 16), CS_dTdX(:, 16), 'k',...
     [0.6, 0.6], [-10, 10], 'r--',...
     [0.4, 0.4], [-10, 10], 'r--',...
     'linewidth', linewdith, 'markersize', markersize)
xlim([0, 1])
ylim([-6, 6])
xlabel('X', 'fontsize', fontsize)
ylabel('dT/dR', 'fontsize', fontsize)
legend('CSA', 'Old CSA Method', 'Complex Step')
% legend('CSA', 'Old CSA Method', 'Complex Step', 'location', 'best')
set(gca, 'fontsize', fontsize)

figure,
plot(X(:, 13), CSA_dTdX(:, 13), 'ro',...
     X(:, 13), CSA_dTdX_Old(:, 13), 'b+',...
     X(:, 13), CS_dTdX(:, 13), 'k',...
     [0.6, 0.6], [-10, 10], 'r--',...
     [0.4, 0.4], [-10, 10], 'r--',...
     'linewidth', linewdith, 'markersize', markersize)
xlim([0, 1])
ylim([-2.5, 2.5])
xlabel('X', 'fontsize', fontsize)
ylabel('dT/dR', 'fontsize', fontsize)
legend('CSA', 'Old CSA Method', 'Complex Step')
% legend('CSA', 'Old CSA Method', 'Complex Step', 'location', 'best')
set(gca, 'fontsize', fontsize)