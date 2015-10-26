clc
clear all
close all
format short g
% ------------------------------------------ %
x = linspace(-2, 2, 20);
y = linspace(-2, 2, 20);
dx = x(2) - x(1);
dy = y(2) - y(1);
[Y, X] = meshgrid(y, x);
% xloc = 0.2 * cos(linspace(0, 2*pi, 100));
% yloc = 0.2 * sin(linspace(0, 2*pi, 100));;
% M = 0;
% Z = X.^2 + Y.^2;
% foo = zeros(length(xloc), 3);
% for iloc = 1:length(xloc)
%     M = delta(X, xloc(iloc), x(2) - x(1)) .* delta(Y, yloc(iloc), y(2) -y(1));
%     foo(iloc, :) = [xloc(iloc), yloc(iloc), trapz(y, trapz(x, M.*Z))];
% end
% % f = phi(X, xloc, x(2) - x(1));
% % trapz(y, trapz(x, M.*Z)) % Double integral
% figure,
% scatter(foo(:,1), foo(:,2), foo(:,3))

M = dDeltadX0(X, 0, dx) .* dDeltadX0(Y, 0, dy);
trapz(y, trapz(x, M))

figure,
surf(X, Y, M)