clc;
clear all;
close all;
format short g;
% ------------------------------------------------------------ %
n = 100;
dt = 0.00001;
x = linspace(0, 1, n+2);
y = linspace(0, 1, n+2);
h = x(2) - x(1);

% xloc = 0.5;
% yloc = 0.5;
% [Y, X] = meshgrid(y(2:end-1), x(2:end-1));
% M = delta(X, xloc, h).*delta(Y, yloc, h);
nloc = 100;
[Y, X] = meshgrid(y(2:end-1), x(2:end-1));
xloc = 0.5 + 0.1 * cos(linspace(0, 2*pi, nloc));
yloc = 0.5 + 0.1 * sin(linspace(0, 2*pi, nloc));
% trapz(y, trapz(x, f.*Z)) % Double integral

Tn = zeros(n*n, 1) + eps;
F = Tn;
t = 0:dt:1;

L = laplacian(n, h);
alpha = -100;
beta = 0;
U = 1;
f_time = 0;
for it=1:length(t)
    RHS = (eye(n*n, n*n) + L * dt * 0.5) * Tn + F;
    Tnp1 = (eye(n*n, n*n) - L * dt * 0.5) \ RHS;
    err = max(max(abs(Tnp1 - Tn) ./ abs(Tn)));
    Tnp1 = reshape(Tnp1,[n,n]);
    Tn = reshape(Tnp1,[n,n]);
    
    unp1 = zeros(nloc, 1);
    un = zeros(nloc, 1);
    for iloc = 1:length(xloc)
        M = delta(X, xloc(iloc), x(2) - x(1)) .* delta(Y, yloc(iloc), y(2) -y(1));
        unp1(iloc) = trapz(y(2:end-1), trapz(x(2:end-1), M.*Tnp1));
        un(iloc) = trapz(y(2:end-1), trapz(x(2:end-1), M.*Tn));
    end
    f_time = f_time + (unp1 - U + un - U) * dt * 0.5;
    f_cont = unp1 - U;
    f = alpha * f_time + beta * f_cont;
	F = 0;
    for iloc = 1:length(xloc)
        M = delta(X, xloc(iloc), x(2) - x(1)) .* delta(Y, yloc(iloc), y(2) -y(1));
        F = F + f(iloc) * M;
    end
%     unp1 = trapz(y(2:end-1), trapz(x(2:end-1), M.*Tnp1));
%     un = trapz(y(2:end-1), trapz(x(2:end-1), M.*Tn));
%     f_time = f_time + (unp1 - U + un - U) * dt * 0.5;
%     f_cont = unp1 - U;
%     f = alpha * f_time + beta * f_cont;
%     f = f * M;
    
%     f = reshape(f, n*n, []);
    F = reshape(F, n*n, []);
    Tnp1 = reshape(Tnp1, n*n, []);
    Tn = reshape(Tn, n*n, []);
    Tn = Tnp1;
    
    T = zeros(n+2, n+2);
    T(2:end-1,2:end-1) = reshape(Tn, n, n);
    [Yplt, Xplt] = meshgrid(y, x);
    figure(1),
    surf(Xplt, Yplt, T)
    zlim([0,1])
    pause(0.1)
    
    if (err < 1e-5) && (it>100)
        break
    else
        disp(['err = ', num2str(err)]);
    end
end

