function [dx] = COM_plan(t,x, cog_h, hstep_leng,x0,domain_type, domains,ndomains)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all
clc

if (nargin<=2)
  % Half step lenth 
  hstep_leng = 0.1;  
  cog_h = 0.6528;
  %Inital condition: COM Position, Velocity and ZMP position/COM
  %acceleration
  x=[0;0;0];
  domain_type = 1;
  Sstime = 2;%domains{1,1}.time;
  Dstime = 1.0;%domains{2,1}.time;
  domains{1,1}.time=2;
  domains{2,1}.time=1.0;
  ndomains=2;
  ST=0.1;
  %% new
  t=0;
end
StepT = 3;
SysTime=0.01;
%% others
xstack=[];
chkstack=[];
nstack=[];
for t=linspace(0,3,3/SysTime);
%% ZMP 
    % hstep_leng=0.1;
    % cog_h=Iset.y0(2);
    % x0=X0;
% Time Specs
% ST = 0.1; % sampling time / time step
    % StepT = 0;
    % for i = 1:ndomains
    %     StepT=StepT+domains{i,1}.time;
    % end
    % if domain.type==1
    currentT=t;
    % else
    % currentT=t+2;
    % end
LfowardT = StepT; % look forward time
StepCount = floor(StepT/ST);
Sstime = domains{1,1}.time;
Dsitme = domains{2,1}.time;
N1 = round(Sstime/StepT*StepCount);
N2 = StepCount-N1;
CycleTime = mod(currentT,StepT);
CycleTimeCount= round(CycleTime/StepT*StepCount);
if CycleTimeCount<1&&currentT<3
    CycleTimeCount=1;
elseif CycleTimeCount<1&&currentT>2.8
    CycleTimeCount=StepCount;  
end
if CycleTimeCount<N1
    ModN1=N1-CycleTimeCount;
    ModN2=N2;
    NxtN1=CycleTimeCount;
    NxtN2=0;
%     StepSelection = blkdiag(ones(N1-CycleTimeCount,1),ones(N2,1),ones(ModN1,1))
elseif CycleTimeCount>=N1;
    ModN1=0;
    NxtN1=N1;
    NxtN2=CycleTimeCount-N1;
    ModN2=N2-NxtN2;
%     StepSelection = blkdiag(ones(N2-NxtN2,1),ones(ModN1,1),ones(NxtN2,1))
end
StepSelection = blkdiag(ones(ModN1,1),ones(N2-NxtN2,1),ones(NxtN1,1),ones(NxtN2,1));
Step1 = StepSelection(:,1);
Step2 = StepSelection(:,2);
Step3 = StepSelection(:,3);
Step4 = StepSelection(:,4);


%% desire terminal com velocity
midV =(hstep_leng/StepT);
xdata = linspace(0,1,StepCount)';
% p2 = FiveOrderTraj(0, 0, 0, 0, 0,0, midV,0.5);
% desireV=polyval(p2,xdata);
% ShiftDesV=[desireV(CycleTimeCount:end);desireV(1:CycleTimeCount-1) ];
ShiftDesV=0.0;
%% desire terminal com pos
desirePos =hstep_leng+CycleTimeCount/StepCount*(hstep_leng);%( pos(2)+pos(12) )/2;
%%
p3 = FiveOrderTraj(0, 0, 0, hstep_leng, 0,0, hstep_leng/2,0.5);
desireZMP=polyval(p3,xdata);
ShiftZMP=[desireZMP(CycleTimeCount:end);desireZMP(1:CycleTimeCount-1)+0.1 ];
%%

%%% LIPM
% Specs
T=ST;
N=StepCount;
g=9.8;
hcom=cog_h;
z0=cog_h;
omega2= g/hcom;
A = [1        T     0;
     omega2*T 1 -omega2*T;
     0        0     1;
     ];
B=[0;0;T];
C=[0 0 1 ]; 
Asys = [1        SysTime     0;
     omega2*SysTime 1 -omega2*SysTime;
     0              0        1;
     ]; 
Bsys=[0;0;SysTime];


AB = [A,B,zeros(3,N-1)];
CC=C*AB;
for i = 1:N-1
    ABtemp = [A*AB(end-2:end,1:3+i),B,zeros(3,N-1-i)];
    AB=[AB;ABtemp];
    CC=[CC;C*ABtemp];
end

Acom = AB(1:3:end,1:3);
Bcom = AB(1:3:end,4:end);
AcomV = AB(2:3:end,1:3);
BcomV = AB(2:3:end,4:end);
Azmp = CC(:,1:3);
Bzmp = CC(:,4:end);

%%% QP
%%
w1=1e-10; %c input
w2=1e-2; % com 
w3=1e-2; % comV
w4=1e0;
%%

% Hzmp =  w1*eye(N)+w2*Bcom'*Bcom+w3*BcomV'*BcomV+w4*Bzmp'*Bzmp;
% %%Insert Terminal COM position (hstep_leng) & Vel (desireV)
% fzmp = w2*Bcom'*(Acom*x-desirePos)+w3*BcomV'*(AcomV*x-ShiftDesV)+w4*Bzmp'*(Azmp*x-ShiftZMP);

Hzmp =  w1*eye(N)+w4*Bzmp'*Bzmp;
%%Insert Terminal COM position (hstep_leng) & Vel (desireV)
fzmp = w4*Bzmp'*(Azmp*x-ShiftZMP);
% ineqaulity constraint
Pzmp1 =  0.15;
Nzmp1 = -0.05;
Pzmp2 =  0.15 + hstep_leng;
Nzmp2 = -0.05;
Pzmp3 =  0.15 + hstep_leng;
Nzmp3 = -0.05 ;
Pzmp4 =  0.15+hstep_leng*2;
Nzmp4 = -0.05+hstep_leng ;
% Ax_ineq = [];
% bx_ineq = [];

Azmp_Piq = Bzmp; bzmp_Piq = -Azmp*x+Pzmp1*Step1+Pzmp2*Step2+Pzmp3*Step3+Pzmp4*Step4;
Azmp_Niq = -Bzmp; bzmp_Niq = Azmp*x-Nzmp1*Step1-Nzmp2*Step2-Nzmp3*Step3-Nzmp4*Step4;

Ax_ineq = [Azmp_Piq;Azmp_Niq];
bx_ineq = [bzmp_Piq;bzmp_Niq];

% Ax_eq = [];
% bx_eq = [];

temp = Acom*x;
Acom_eq = Bcom(end,:); bcom_eq = desirePos-temp(end);
temp = AcomV*x;
AcomV_eq = BcomV(end,:); bcomV_eq = ShiftDesV(end)-temp(end);

ComBound=0.01;
ComVBound=0.01;
temp = Acom*x;
Acom_Piq = Bcom(end,:); bcom_Piq = desirePos+ComBound-temp(end);
Acom_Niq = -1*Bcom(end,:); bcom_Niq = temp(end)-(desirePos-ComBound);
temp = AcomV*x;
AcomV_Piq = BcomV(end,:); bcomV_Piq = ShiftDesV(end)+ComVBound-temp(end);
AcomV_Niq = -1*BcomV(end,:); bcomV_Niq = temp(end)-(ShiftDesV(end)-ComVBound);

Ax_ineq=[Ax_ineq;Acom_Piq;Acom_Niq;AcomV_Piq;AcomV_Niq];
bx_ineq=[bx_ineq;bcom_Piq;bcom_Niq;bcomV_Piq;bcomV_Niq];


Ax_eq=[Acom_eq;AcomV_eq];
bx_eq=[bcom_eq;bcomV_eq];
% Ax_eq = zeros(size([AcomV_eq]));
% bx_eq = zeros(size([bcomV_eq]));
 
options = optimset('Algorithm','interior-point-convex','Display','off');
% 
u0 = -Hzmp\fzmp;
% 
% save('test.mat','H','f')
% 
Hzmp=(Hzmp+Hzmp')/2;
[u1,fval,existflag] = quadprog(Hzmp,fzmp,[],[],Ax_eq,bx_eq,[],[],u0,options);%Ax_ineq bx_ineq
    fprintf('\ntime = %e, existflag = %u , fval = %e',t,existflag,fval)  ;     %existflag

    Acon = [0 1 0 ;
         omega2 0 -omega2 ;
         0 0 0];
    Bcon = [0;0;1] ;
    dx=Acon*x+Bcon*u1(1)/10^6;
% Z = Azmp*x+Bzmp*u1;
% C = Acom*x+Bcom*u1;
% CV = AcomV*x+BcomV*u1;
    noise=0;randn(3,1)*10^-3;
    nstack=[nstack noise];
    x=Asys*x+Bsys*u1(1)+noise;
    xstack=[xstack x];
    chkstack=[chkstack ShiftZMP];
% % % 
% %% Curve fitting
% xdata = linspace(0,1,StepCount)';
% ydata = C;
% PolyOrder=7;
% p = mmpolyfit(xdata,ydata,PolyOrder,'Point',[xdata(1) ydata(1);xdata(end) ydata(end)],'Slope',[xdata(1:2) zeros(2,1);xdata(end-1:end) zeros(2,1)]);
% ydataLearned = polyval(p,xdata);
% 
% % p2 = FiveOrderTraj(0, 0, 0, 0.001, 0,hstep_leng/StepT, 0.0);
% % xdata2 = linspace(0,1,StepCount*10)';
% % dataComV=polyval(p2,xdata);
% h=figure();
% % set(h,'position',[100,100,1024*1.5,768*1.5]);
% hold on
% % plot(Zmprefx(1:N),'r','Linewidth',3,'LineStyle','-.')
% plot(Z,'Linewidth',2)
% plot(C,'g','Linewidth',2)
% plot(ydataLearned,'m','Linewidth',2)
% plot(CV,'kx','Markersize',3)
% plot(ShiftDesV,'g','Linewidth',5)
% hold off
% legend('ZMP from QP','COG from QP', 'Polyfit-7 order','COG Velocity')



end
%% Save as a function

% syms x t real
% % A = sym('A', [3 4])
% for i=1:PolyOrder+1
%     if i==1
%     x = p(PolyOrder+2-i);
%     tPoly = t/StepT;
%     else
%     x = x+p(PolyOrder+2-i)*tPoly;
%     tPoly = tPoly*t/StepT;
%     end
%     
% end
% 
% x = simplify(x);
% dxdt = simplify(diff(x,t));
% ddxdt = simplify(diff(dxdt,t));
% 
% matlabFunction(x, 'file', 'ZMPTrajGenX_v2');
% matlabFunction(dxdt, 'file', 'ZMPTrajGenXder_v2');
% matlabFunction(ddxdt, 'file', 'ZMPTrajGenXder2_v2');

% PolyFuncGen(PolyOrder,StepT)
xdata = linspace(0,1,StepCount*ST/SysTime)';
desireZMP=polyval(p3,xdata);
plot(xstack','DisplayName','xstack') 
hold on
plot(desireZMP,'g','Linewidth',3)
hold off

% plot(chkstack)
end

function [p]= FiveOrderTraj(startP,startV,startA,endP,endV,endA,midP,midT)
% p(6)=startP;
% p(5)=startV;
% p(4)=startA/2;
% p(3)=(20*endP-20*startP-(8*endV+12*startV)*endT-(3*startA-endA)*endT^2)/2/endT^3;
% P(2)=(30*startP-30*endP+(14*endV+16*startV)*endT+(3*startA-2*endA)*endT^2)/2/endT^4;
% p(1)=(12*endP-12*startP-(6*endV+6*startV)*endT-(startA-endA)*endT^2 )/2/endT^5;
% xstart=0;
PolyOrder = 6;
p=zeros(PolyOrder+1,1);
T=[ 1 0 0 0 0 0 0;
    0 1 0 0 0 0 0;
    0 0 2 0 0 0 0;
    1 1 1 1 1 1 1;
    0 1 2 3 4 5 6;
    0 0 2 6 12 20 30;    
    1 midT midT^2 midT^3 midT^4 midT^5 midT^6;
  ];
A=T^-1*[startP;startV;startA;endP;endV;endA;midP];
for i=1:PolyOrder+1
    p(i)=A(PolyOrder+1-i+1);
end

end
