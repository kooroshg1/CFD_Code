clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
n = 50;
dt = 1e-2;
x = linspace(0, 1, n+2);
y = linspace(0, 1, n+2);
h = x(2) - x(1);

xloc = 0.5 + 0.001i;
yloc = 0.5;
[Y, X] = meshgrid(y(2:end-1), x(2:end-1));
M = delta(X, xloc, h).*delta(Y, yloc, h);
% trapz(y, trapz(x, f.*Z)) % Double integral
% break
Tn = sparse(n*n, 1) + eps;
f = Tn;
t = 0:dt:1;

L = laplacian(n, h);
alpha = -100;
beta = -0;
U = 1;
f_time = 0;
for it=1:length(t)
    RHS = (speye(n*n, n*n) + L * dt * 0.5) * Tn + f * dt;
    Tnp1 = (speye(n*n, n*n) - L * dt * 0.5) \ RHS;
    
    err = max(max(abs(Tnp1 - Tn) ./ abs(Tn)));
    Tnp1 = reshape(Tnp1,[n,n]);
    Tn = reshape(Tn,[n,n]);

    unp1 = trapz(y(2:end-1), trapz(x(2:end-1), M.*Tnp1));
    un = trapz(y(2:end-1), trapz(x(2:end-1), M.*Tn));
    f_time = f_time + (unp1 - U + un - U) * dt * 0.5;
    f_cont = unp1 - U;
    f = alpha * f_time + beta * f_cont;
    f = f * M;
    
    f = reshape(f, n*n, []);
    Tnp1 = reshape(Tnp1, n*n, []);
    Tn = Tnp1;
    
    T = zeros(n+2, n+2);
    T(2:end-1,2:end-1) = reshape(Tn, n, n);
    [Yplt, Xplt] = meshgrid(y, x);
%     figure(1),
%     surf(Xplt, Yplt, T)
%     zlim([0,1])
    
    if (err < 1e-5) && (it>100)
        break
    elseif rem(it,10) == 0
        disp(['err = ', num2str(err)]);
    end
end
dlmwrite('X.txt', Xplt)
dlmwrite('Y.txt', Yplt)
dlmwrite('T1.txt',T)
