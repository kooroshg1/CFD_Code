clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
xImm100X100 = dlmread('100X100/xImm');
yImm100X100 = dlmread('100X100/yImm');
Xu100X100 = dlmread('100X100/Xu.txt');
Yu100X100 = dlmread('100X100/Yu.txt');
Xp100X100 = dlmread('100X100/Xp.txt');
Yp100X100 = dlmread('100X100/Yp.txt');
U100X100_CSA = dlmread('100X100/Up.txt');
V100X100_CSA = dlmread('100X100/Vp.txt');
P100X100_CSA = dlmread('100X100/Pp.txt');
pBoundary100X100_CSA = dlmread('100X100/ppBoundary.txt');

U100X100_CS = dlmread('100X100/U_CS.txt');
V100X100_CS = dlmread('100X100/V_CS.txt');
P100X100_CS = dlmread('100X100/P_CS.txt');
pBoundary100X100_CS = dlmread('100X100/pBoundary_CS.txt');

[Uin100,Uon100] = inpolygon(Xu100X100, Yu100X100, xImm100X100, yImm100X100);
U100X100_CS(Uin100) = nan;
U100X100_CSA(Uin100) = nan;
% ----------------------------------------------------------------------- %
xImm200X200 = dlmread('200X200/xImm');
yImm200X200 = dlmread('200X200/yImm');
Xu200X200 = dlmread('200X200/Xu.txt');
Yu200X200 = dlmread('200X200/Yu.txt');
Xp200X200 = dlmread('200X200/Xp.txt');
Yp200X200 = dlmread('200X200/Yp.txt');
U200X200_CSA = dlmread('200X200/Up.txt');
V200X200_CSA = dlmread('200X200/Vp.txt');
P200X200_CSA = dlmread('200X200/Pp.txt');
pBoundary200X200_CSA = dlmread('200X200/ppBoundary.txt');

U200X200_CS = dlmread('200X200/U_CS.txt');
V200X200_CS = dlmread('200X200/V_CS.txt');
P200X200_CS = dlmread('200X200/P_CS.txt');
pBoundary200X200_CS = dlmread('200X200/pBoundary_CS.txt');

[Uin200,Uon200] = inpolygon(Xu200X200, Yu200X200, xImm200X200, yImm200X200);
U200X200_CS(Uin200) = nan;
U200X200_CSA(Uin200) = nan;
% ----------------------------------------------------------------------- %
xImm300X300 = dlmread('300X300/xImm');
yImm300X300 = dlmread('300X300/yImm');
Xu300X300 = dlmread('300X300/Xu.txt');
Yu300X300 = dlmread('300X300/Yu.txt');
Xp300X300 = dlmread('300X300/Xp.txt');
Yp300X300 = dlmread('300X300/Yp.txt');
U300X300_CSA = dlmread('300X300/Up.txt');
V300X300_CSA = dlmread('300X300/Vp.txt');
P300X300_CSA = dlmread('300X300/Pp.txt');
pBoundary300X300_CSA = dlmread('300X300/ppBoundary.txt');

U300X300_CS = dlmread('300X300/U_CS.txt');
V300X300_CS = dlmread('300X300/V_CS.txt');
P300X300_CS = dlmread('300X300/P_CS.txt');
pBoundary300X300_CS = dlmread('300X300/pBoundary_CS.txt');

[Uin300,Uon300] = inpolygon(Xu300X300, Yu300X300, xImm300X300, yImm300X300);
U300X300_CS(Uin300) = nan;
U300X300_CSA(Uin300) = nan;
% ----------------------------------------------------------------------- %
xImm400X400 = dlmread('400X400/xImm');
Xu400X400 = dlmread('400X400/Xu.txt');
Yu400X400 = dlmread('400X400/Yu.txt');
Xp400X400 = dlmread('400X400/Xp.txt');
Yp400X400 = dlmread('400X400/Yp.txt');
U400X400_CSA = dlmread('400X400/Up.txt');
V400X400_CSA = dlmread('400X400/Vp.txt');
P400X400_CSA = dlmread('400X400/Pp.txt');
pBoundary400X400_CSA = dlmread('400X400/ppBoundary.txt');

U400X400_CS = dlmread('400X400/U_CS.txt');
V400X400_CS = dlmread('400X400/V_CS.txt');
P400X400_CS = dlmread('400X400/P_CS.txt');
pBoundary400X400_CS = dlmread('400X400/pBoundary_CS.txt');


yint = 0.65;
[val, ind100] = min(abs(Xp100X100(:,1)' - yint));
[val, ind200] = min(abs(Xp200X200(:,1)' - yint));
[val, ind300] = min(abs(Xp300X300(:,1)' - yint));
[val, ind400] = min(abs(Xp400X400(:,1)' - yint));

figure,
plot(Xp200X200(:,ind200), P200X200_CS(:,ind200),'r',...
     Xp200X200(:,ind200), P200X200_CSA(:,ind200),'k')
 
figure,
plot(Xp300X300(:,ind300), P300X300_CS(:,ind300),'r',...
     Xp300X300(:,ind300), P300X300_CSA(:,ind300),'k')
 
figure,
plot(Xp400X400(:,ind400), P400X400_CS(:,ind400),'r',...
     Xp400X400(:,ind400), P400X400_CSA(:,ind400),'k')

% figure,
% plot(Yu300X300(ind300,:), U300X300_CS(ind300,:),'r',...
%      Yu300X300(ind300,:), U300X300_CSA(ind300,:),'k')
% [val, ind100] = min(abs(Xu100X100(:,1)' - yint));
% [val, ind200] = min(abs(Xu200X200(:,1)' - yint));
% [val, ind300] = min(abs(Xu300X300(:,1)' - yint));
% [val, ind400] = min(abs(Xu400X400(:,1)' - yint));
% figure,
% plot(Xu100X100(:,ind100), U100X100_CSA(:,ind100),...
%      Xu200X200(:,ind200), U200X200_CSA(:,ind200),...
%      Xu300X300(:,ind300), U300X300_CSA(:,ind300),...
%      Xu400X400(:,ind400), U400X400_CSA(:,ind400))
% legend('100','200','300','400') 
% 
% figure,
% plot(Xu100X100(:,ind100), U100X100_CS(:,ind100),...
%      Xu200X200(:,ind200), U200X200_CS(:,ind200),...
%      Xu300X300(:,ind300), U300X300_CS(:,ind300),...
%      Xu400X400(:,ind400), U400X400_CS(:,ind400))
% legend('100','200','300','400') 
% 

xint = 0.65;
[val, ind100] = min(abs(Yu100X100(1,:)' - xint));
[val, ind200] = min(abs(Yu200X200(1,:)' - xint));
[val, ind300] = min(abs(Yu300X300(1,:)' - xint));
[val, ind400] = min(abs(Yu400X400(1,:)' - xint));

figure,
plot(Yu200X200(ind200,:), U200X200_CS(ind200,:),'r',...
     Yu200X200(ind200,:), U200X200_CSA(ind200,:),'k')

figure,
plot(Yu300X300(ind300,:), U300X300_CS(ind300,:),'r',...
     Yu300X300(ind300,:), U300X300_CSA(ind300,:),'k')

% figure,
% plot(xImm100X100, pBoundary100X100_CS,'k--',...
%      xImm200X200, pBoundary200X200_CS,'r--',...
%      xImm300X300, pBoundary300X300_CS,'b--')
% hold on
% plot(xImm100X100, pBoundary100X100_CSA,'k',...
%      xImm200X200, pBoundary200X200_CSA,'r',...
%      xImm300X300, pBoundary300X300_CSA,'b')
% legend('100','200','300') 