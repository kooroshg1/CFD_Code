clc;
clear all;
close all;
format short g;
% ------------------------------------------------------------ %
n = 30;
dt = 1e-3;
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

% ANALYSIS
Tn = sparse(n*n, 1) + eps;
F = Tn;
f_time = 0;
U = 1;
% SA
Tsn = sparse(n*n, 1) + eps;
Fs = Tsn;
fs_time = 0;

t = 0:dt:1;

L = laplacian(n, h);
alpha = -100;
beta = -0;
for it=1:length(t)
    % ANALYSIS
    RHS = (speye(n*n, n*n) + L * dt * 0.5) * Tn + F * dt;
    Tnp1 = (speye(n*n, n*n) - L * dt * 0.5) \ RHS;
    % SA
    RHSs = (speye(n*n, n*n) + L * dt * 0.5) * Tsn + Fs * dt;
    Tsnp1 = (speye(n*n, n*n) - L * dt * 0.5) \ RHSs;
    
    % ANALYSIS
    err = max(max(abs(Tnp1 - Tn) ./ abs(Tn)));
    % SA
    errs = max(max(abs(Tsnp1 - Tsn) ./ abs(Tsn)));
    
    % ANALYSIS
    Tnp1 = reshape(Tnp1,[n,n]);
    Tn = reshape(Tn,[n,n]);
    % SA
    Tsnp1 = reshape(Tsnp1,[n,n]);
    Tsn = reshape(Tsn,[n,n]);
    
    % ANALYSIS
    unp1 = zeros(nloc, 1);
    un = zeros(nloc, 1);
    % SA
    usnp1 = zeros(nloc, 1);
    usn = zeros(nloc, 1);
    
    
    for iloc = 1:length(xloc)
        % ANALYSIS
        M = delta(X, xloc(iloc), x(2) - x(1)) .* delta(Y, yloc(iloc), y(2) -y(1));
        unp1(iloc) = trapz(y(2:end-1), trapz(x(2:end-1), M.*Tnp1));
        un(iloc) = trapz(y(2:end-1), trapz(x(2:end-1), M.*Tn));
        % SA
        dMdX = dDeltadX0(X, xloc(iloc), x(2) - x(1)) .* delta(Y, yloc(iloc), y(2) -y(1));
        usnp1(iloc) = trapz(y(2:end-1), trapz(x(2:end-1), M .* Tsnp1));
        usnp1(iloc) = usnp1(iloc) + trapz(y(2:end-1), trapz(x(2:end-1), dMdX .* Tnp1));
        
        usn(iloc) = trapz(y(2:end-1), trapz(x(2:end-1), M .* Tsn));
        usn(iloc) = usn(iloc) + trapz(y(2:end-1), trapz(x(2:end-1), dMdX .* Tn));
    end

    % ANALYSIS
    f_time = f_time + (unp1 - U + un - U) * dt * 0.5;
    f_cont = unp1 - U;
    f = alpha * f_time + beta * f_cont;
	F = 0;
    % SA
    fs_time = fs_time + (usnp1 + usn) * dt * 0.5;
    fs_cont = usnp1;
    fs = alpha * fs_time + beta * fs_cont;
	Fs = 0;

    for iloc = 1:length(xloc)
        % ANALYSIS
        M = delta(X, xloc(iloc), x(2) - x(1)) .* delta(Y, yloc(iloc), y(2) -y(1));
        F = F + f(iloc) * M;
        % SA
        dMdX = dDeltadX0(X, xloc(iloc), x(2) - x(1)) .* delta(Y, yloc(iloc), y(2) -y(1));
        Fs = Fs + fs(iloc) * M + f(iloc) * dMdX;
    end
    
    % ANALYSIS
    F = reshape(F, n*n, []);
    Tnp1 = reshape(Tnp1, n*n, []);
    Tn = Tnp1;
    % SA
    Fs = reshape(Fs, n*n, []);
    Tsnp1 = reshape(Tsnp1, n*n, []);
    Tsn = Tsnp1;
    
    % ANALYSIS
    T = zeros(n+2, n+2);
    T(2:end-1,2:end-1) = reshape(Tn, n, n);
    % SA
    Ts = zeros(n+2, n+2);
    Ts(2:end-1,2:end-1) = reshape(Tsn, n, n);
    [Yplt, Xplt] = meshgrid(y, x);
    figure(1),
    surf(Xplt, Yplt, T, 'linestyle', 'none')
    zlim([0,1])
    figure(2),
    surf(Xplt, Yplt, Ts, 'linestyle', 'none')
%     zlim([0,1])
    if (err < 1e-5) && (it>100)
        break
    elseif rem(it,10) == 0
        disp(['err = ', num2str(err)]);
    end
end
dlmwrite('X.txt', Xplt)
dlmwrite('Y.txt', Yplt)
dlmwrite('T0.txt',T)
dlmwrite('CSA_dTdX.txt',Ts)

