clc;
% clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
fontsize = 40;
linewidth = 4;
% ----------------------------------------------------------------------- %
xImm = dlmread('xImm');
yImm = dlmread('yImm');
Xu = dlmread('Xu.txt');
Yu = dlmread('Yu.txt');
Xv = dlmread('Xv.txt');
Yv = dlmread('Yv.txt');
Xp = dlmread('Xp.txt');
Yp = dlmread('Yp.txt');

% pBoundary = dlmread('plongBoundary.txt');
% U = dlmread('Ulong.txt');
% V = dlmread('Vlong.txt');
% P = dlmread('Plong.txt');

p0Boundary = dlmread('p0Boundary.txt');
U0 = dlmread('U0.txt');
V0 = dlmread('V0.txt');
P0 = dlmread('P0.txt');

p1Boundary = dlmread('p1Boundary.txt');
U1 = dlmread('U1.txt');
V1 = dlmread('V1.txt');
P1 = dlmread('P1.txt');

p2Boundary = dlmread('p2Boundary.txt');
U2 = dlmread('U2.txt');
V2 = dlmread('V2.txt');
P2 = dlmread('P2.txt');

% ppBoundary = dlmread('ppBoundary.txt');
% Up = dlmread('Up.txt');
% Vp = dlmread('Vp.txt');
% Pp = dlmread('Pp.txt');

[Uin, Uon] = inpolygon(Xu, Yu, xImm, yImm); 
[Vin, Von] = inpolygon(Xv, Yv, xImm, yImm);
[Pin, Pon] = inpolygon(Xp, Yp, xImm, yImm);

U0(Uin) = nan; U1(Uin) = nan; U0(Uon) = nan; U1(Uon) = nan;
V0(Vin) = nan; V1(Vin) = nan; V0(Von) = nan; V1(Von) = nan;
P0(Pin) = nan; P1(Pin) = nan; P0(Pon) = nan; P1(Pon) = nan;
P(Pin) = nan;
U(Uin) = nan;
V(Vin) = nan;

step_size1 = 0.001;
step_size2 = 0.0001;
Pp_FD1 = (P1 - P0) / step_size1;
Up_FD1 = (U1 - U0) / step_size1;
Vp_FD1 = (V1 - V0) / step_size1;
ppBoundary_FD1 = (p1Boundary - p0Boundary) / step_size1;

Pp_FD2 = (P2 - P0) / step_size2;
Up_FD2 = (U2 - U0) / step_size2;
Vp_FD2 = (V2 - V0) / step_size2;
ppBoundary_FD2 = (p2Boundary - p0Boundary) / step_size2;

% % Up(Uin) = nan;
% % Vp(Vin) = nan;
% % Pp(Pin) = nan;

figure,
contourf(Xu, Yu, U, 50, 'linestyle', 'none')
colorbar('northoutside')
xlabel('X','fontsize',fontsize)
ylabel('Y','fontsize',fontsize)
set(gca,'fontsize',fontsize)

figure,
contourf(Xp, Yp, P, 50, 'linestyle', 'none')
colorbar('northoutside')
xlabel('X','fontsize',fontsize)
ylabel('Y','fontsize',fontsize)
set(gca,'fontsize',fontsize)

figure,
contourf(Xu, Yu, Up_FD1, 50, 'linestyle', 'none')
colorbar('northoutside')
xlabel('X','fontsize',fontsize)
ylabel('Y','fontsize',fontsize)
set(gca,'fontsize',fontsize)

figure,
contourf(Xu, Yu, Up_FD2, 50, 'linestyle', 'none')
colorbar('northoutside')
xlabel('X','fontsize',fontsize)
ylabel('Y','fontsize',fontsize)
set(gca,'fontsize',fontsize)

% figure,
% contourf(Xp, Yp, Pp_FD1, 50, 'linestyle', 'none')
% 
% figure,
% contourf(Xp, Yp, Pp_FD2, 50, 'linestyle', 'none')


filterDegree = 7;
filterSize = 41;
figure,
plot(xImm, sgolayfilt(ppBoundary_FD1,filterDegree,filterSize),'r',...
     xImm, sgolayfilt(ppBoundary_FD2,filterDegree,filterSize),'k',...
     'linewidth', linewidth)
xlabel('Chord', 'fontsize', fontsize)
ylabel('dP/db', 'fontsize', fontsize)
legend('CSA', 'CS')
set(gca,'fontsize',fontsize)
 
trapz(xImm, ppBoundary_FD1)
trapz(xImm, ppBoundary_FD2)