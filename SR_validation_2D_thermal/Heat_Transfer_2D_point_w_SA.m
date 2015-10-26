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

xloc = 0.5;
yloc = 0.5;
[Y, X] = meshgrid(y(2:end-1), x(2:end-1));
M = delta(X, xloc, h) .* delta(Y, yloc, h);
dMdX = dDeltadX0(X, xloc, h) .* delta(Y, yloc, h);
dMdY = delta(X, xloc, h) .* dDeltadX0(Y, yloc, h);
% trapz(y, trapz(x, f.*Z)) % Double integral

% ANALYSIS
Tn = sparse(n*n, 1) + eps;
f = Tn;
U = 1;
f_time = 0;
% SA
Tsn = sparse(n*n, 1) + eps;
fs = Tsn;
fs_time = 0;

t = 0:dt:1;

L = laplacian(n, h);
alpha = -100;
beta = -0;
for it=1:length(t)
    % ANALYSIS
    RHS = (speye(n*n, n*n) + L * dt * 0.5) * Tn + f * dt;
    Tnp1 = (speye(n*n, n*n) - L * dt * 0.5) \ RHS;
    % SA
    RHSs = (speye(n*n, n*n) + L * dt * 0.5) * Tsn + fs * dt;
    Tsnp1 = (speye(n*n, n*n) - L * dt * 0.5) \ RHSs;
    
    % ANALYSIS
    err = max(max(abs(Tnp1 - Tn) ./ abs(Tn)));
    Tnp1 = reshape(Tnp1,[n,n]);
    Tn = reshape(Tn,[n,n]);
    % SA
    errs = max(max(abs(Tsnp1 - Tsn) ./ abs(Tsn)));
    Tsnp1 = reshape(Tsnp1,[n,n]);
    Tsn = reshape(Tsnp1,[n,n]);
    
    % ANALYSIS
    unp1 = trapz(y(2:end-1), trapz(x(2:end-1), M.*Tnp1));
    un = trapz(y(2:end-1), trapz(x(2:end-1), M.*Tn));
    f_time = f_time + (unp1 - U + un - U) * dt * 0.5;
    f_cont = unp1 - U;
    f = alpha * f_time + beta * f_cont;
    f = f * M;
    
    % SA
    usnp1 = trapz(y(2:end-1), trapz(x(2:end-1), M .* Tsnp1));
    usnp1 = usnp1 + trapz(y(2:end-1), trapz(x(2:end-1), dMdX .* Tnp1));
    usn = trapz(y(2:end-1), trapz(x(2:end-1), M .* Tsn));
    usn = usn + trapz(y(2:end-1), trapz(x(2:end-1), dMdX .* Tn));
    fs_time = fs_time + (usnp1 + usn) * dt * 0.5;
    fs_cont = unp1 - U;
    fs = alpha * fs_time + beta * fs_cont;
    fs = fs * M + (alpha * f_time + beta * f_cont) * dMdX;
    
    % ANALYSIS
    f = reshape(f, n*n, []);
    Tnp1 = reshape(Tnp1, n*n, []);
    Tn = Tnp1;
    % SA
    fs = reshape(fs, n*n, []);
    Tsnp1 = reshape(Tsnp1, n*n, []);
    Tsn = Tsnp1;
    
    % ANALYSIS
    T = zeros(n+2, n+2);
    T(2:end-1,2:end-1) = reshape(Tn, n, n);
    [Yplt, Xplt] = meshgrid(y, x);
    figure(1),
    surf(Xplt, Yplt, T)
    zlim([0,1])
    
    % SA
    Ts = zeros(n+2, n+2);
    Ts(2:end-1,2:end-1) = reshape(Tsn, n, n);
    figure(2),
    surf(Xplt, Yplt, Ts)
%     zlim([0,1])
    
    if (err < 1e-5) && (it>100)
        break
    elseif rem(it,10) == 0
        disp(['err = ', num2str(err)]);
    end
end
% dlmwrite('SA_X.txt', Xplt)
% dlmwrite('SA_Y.txt', Yplt)
% dlmwrite('SA_dTdX.txt',Ts)
