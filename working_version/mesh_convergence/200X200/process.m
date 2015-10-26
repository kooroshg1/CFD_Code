clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
xImm = dlmread('xImm');
yImm = dlmread('yImm');
Xu = dlmread('Xu.txt');
Yu = dlmread('Yu.txt');
Xv = dlmread('Xv.txt');
Yv = dlmread('Yv.txt');
Xp = dlmread('Xp.txt');
Yp = dlmread('Yp.txt');

U0 = dlmread('U0.txt');
V0 = dlmread('V0.txt');
P0 = dlmread('P0.txt');
p0Boundary = dlmread('p0Boundary.txt');

U1 = dlmread('U1.txt');
V1 = dlmread('V1.txt');
P1 = dlmread('P1.txt');
p1Boundary = dlmread('p1Boundary.txt');

U_CSA = dlmread('Up.txt');
V_CSA = dlmread('Vp.txt');
P_CSA = dlmread('Pp.txt');
pBoundary_CSA = dlmread('ppBoundary.txt');

U_CS = dlmread('U_CS.txt');
V_CS = dlmread('V_CS.txt');
P_CS = dlmread('P_CS.txt');
pBoundary_CS = dlmread('pBoundary_CS.txt');

[Uin, Uon] = inpolygon(Xu, Yu, xImm, yImm);
[Vin, Von] = inpolygon(Xv, Yv, xImm, yImm);
[Pin, Pon] = inpolygon(Xp, Yp, xImm, yImm);

U0(Uin) = nan; U1(Uin) = nan;
V0(Vin) = nan; V1(Vin) = nan;
P0(Pin) = nan; P1(Pin) = nan;
U_CSA(Uin) = nan;
V_CSA(Vin) = nan;
P_CSA(Pin) = nan;
U_CS(Uin) = nan;
V_CS(Vin) = nan;
P_CS(Pin) = nan;

U_FD = (U1 - U0) / -0.0001;
V_FD = (V1 - V0) / -0.0001;
P_FD = (P1 - P0) / -0.0001;
pBoundary_FD = (p1Boundary - p0Boundary) / -0.0001;

figure,
subplot(1,2,1)
contourf(Xu,Yu,U_CS,50,'linestyle','none')
axis equal
subplot(1,2,2)
contourf(Xu,Yu,U_CSA,50,'linestyle','none')
axis equal

% figure,
% contourf(Xu,Yu,U_CS,50,'linestyle','none')

yint = 0.2;
[val, ind] = min(abs(Yp(1,:)' - yint))
figure,
plot(Yp(ind,:), P_CS(:,ind),...
     Yp(ind,:), P_CSA(:,ind))
legend('FD','CSA') 


xint = 0.2;
[val, ind] = min(abs(Xu(:,1)' - xint))
figure,
plot(Xu(:,ind), U_CS(:,ind),...
     Xu(:,ind), U_CSA(:,ind))
legend('FD','CSA') 