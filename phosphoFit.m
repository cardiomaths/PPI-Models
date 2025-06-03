function phosphoFit()

%Randomises the LatinHypercube 
rng('shuffle');

% Import any data
dPI = csvread("data/ExperimentalData.csv",2,0);
dataPI = dPI(1:8,2);

%Set time - I need the time points collected to match data time points
%---------------------------------------------------------------------
t=[0,30,60,90,120,180,600];
noTimepoints=7; 
tspan=t; %we want x at every t

%Initial Conditions
%------------------
u1=2700000; %PI
u2=640000;  %PI4
u3=310000;  %PIP2 
u4=1900;    %PIP3 
u5=5200;    %PI34
u6=25000;   %Pp
u7=0; %IP3
u8=0;   %IP1
u9=100000000; %Ipool

u0=[u1 u2 u3 u4 u5 u6 u7 u8 u9];

%Parameters - not needed but included as a reminder of parameter names
par0(1) = 1;  %r1
par0(2) = 1;  %r1m
par0(3) = 1;  %r2
par0(4) = 1;  %r2m
par0(5) = 0.04;  %r3 – we know this parameter
par0(6) = 1;  %r4
par0(7) = 0.0002;  %t1 – we know this parameter
par0(8) = 1;  %t2
par0(9) = 1; %t3
par0(10) = 1; %t3m 
par0(11) = 1; %t4 
par0(12) = 1; %t4m
par0(13) = 1; %t5 
par0(14) = 1; %t5m
par0(15) = 1; %s1
par0(16) = 1; %s2
par0(17) = 1; %s2m

%Set lower and upper bounds for each parameter (fmincon boundaries)
%-------------------------------------------------------------------
lb1=0.001; ub1=100; %r1
lb2=0.001; ub2=100; %r1m
lb3=0.001; ub3=100; %r2
lb4=0.001; ub4=100; %r2
lb5=0.0001; ub5=0.001; %r3 – we know this parameter
lb6=0.001; ub6=100; %r4
lb7=0.01; ub7=0.1; %t1 – we know this parameter
lb8=0.001; ub8=100; %t2 - move out of Ipool to membrane
lb9=0.001; ub9=100; %t3
lb10=0.001; ub10=100; %t3m
lb11=0.0001; ub11=10; %t4 - this parameter is critical at preventing initial sharp peaks (it must be able to go low)
lb12=0.001; ub12=100; %t4m
lb13=0.001; ub13=100; %t5
lb14=0.001; ub14=100; %t5m
lb15=0.0001; ub15=10; %s1 %allowed to go lower than s2
lb16=0.001; ub16=100; %s2
lb17=0.001; ub17=100; %s2m

% Initialise arrays for fmincon and choose options
A=[]; b=[]; Aeq=[]; beq=[]; nonlcon=[];

lb=[lb1 lb2 lb3 lb4 lb5 lb6 lb7 lb8 lb9 lb10 lb11 lb12 lb13 lb14 lb15 lb16 lb17]; %lower bounds for parameters
ub=[ub1 ub2 ub3 ub4 ub5 ub6 ub7 ub8 ub9 ub10 ub11 ub12 ub13 ub14 ub15 ub16 ub17]; %Upper bounds for parameters

%Fmincon options: TolX - stop when diff is less than this value
MFE=1000;       % Max Number of iterations
TX=1e-16;       % Tolerance
options=optimset('Algorithm','interior-point', 'MaxFunEvals', MFE, 'TolX',TX);

fval = zeros(length(t)); exitflag = 0; cnt=0;

%Set up the LatinHypercube, sampling parameter values from this
%--------------------------------------------------------------
LHSlb = [-3 -3 -3 -3 -4 -3 -2 -3 -3 -3 -5 -3 -3 -3 -4 -3 -3]; % lower bounds of parameters 
LHSub = [2 2 2 2 -3 2 -1 2 2 2 0 2 2 2 1 2 2]; % upper bounds of parameters 

nParams = 17;   % Number of parameters/variables (must equal number in lb and ub)
nDivs = 10000;   % Number of divisions of parameter space (larger is better, might be 10000?)
X = lhsdesign(nDivs,nParams,'smooth','off');
%save('LHS1000s1.mat','X');

[iMax,p]=size(X);
%iMax=1;

r1space=logspace(LHSlb(1),LHSub(1),1000);
r1mspace=logspace(LHSlb(2),LHSub(2),1000);
r2space=logspace(LHSlb(3),LHSub(3),1000);
r2mspace=logspace(LHSlb(4),LHSub(4),1000);
r3space=logspace(LHSlb(5),LHSub(5),1000);
r4space=logspace(LHSlb(6),LHSub(6),1000);
t1space=logspace(LHSlb(7),LHSub(7),1000);
t2space=logspace(LHSlb(8),LHSub(8),1000);
t3space=logspace(LHSlb(9),LHSub(9),1000);
t3mspace=logspace(LHSlb(10),LHSub(10),1000);
t4space=logspace(LHSlb(11),LHSub(11),1000);
t4mspace=logspace(LHSlb(12),LHSub(12),1000);
t5space=logspace(LHSlb(13),LHSub(13),1000);
t5mspace=logspace(LHSlb(14),LHSub(14),1000);
s1space=logspace(LHSlb(15),LHSub(15),1000);
s2space=logspace(LHSlb(16),LHSub(16),1000);
s2mspace=logspace(LHSlb(17),LHSub(17),1000);

AllParam = [];
%csvwrite('AllParam.csv',AllParam);

for i=1:iMax
     tic
     par0(1)=r1space(ceil(X(i,1)*1000));
     par0(2)=r1mspace(ceil(X(i,2)*1000));
     par0(3)=r2space(ceil(X(i,3)*1000));
     par0(4)=r2mspace(ceil(X(i,4)*1000));
     par0(5)=r3space(ceil(X(i,5)*1000));
     par0(6)=r4space(ceil(X(i,6)*1000));
     par0(7)=t1space(ceil(X(i,7)*1000));
     par0(8)=t2space(ceil(X(i,8)*1000));
     par0(9)=t3space(ceil(X(i,9)*1000));
     par0(10)=t3mspace(ceil(X(i,10)*1000));
     par0(11)=t4space(ceil(X(i,11)*1000));
     par0(12)=t4mspace(ceil(X(i,12)*1000));
     par0(13)=t5space(ceil(X(i,13)*1000));
     par0(14)=t5mspace(ceil(X(i,14)*1000));
     par0(15)=s1space(ceil(X(i,15)*1000));
     par0(16)=s2space(ceil(X(i,16)*1000));
     par0(17)=s2mspace(ceil(X(i,17)*1000));

     dlmwrite('PriorA0.csv',par0,'-append');
     index=i
     [opt_param,fval,exitflag] = fmincon(@FminconFun,par0,A,b,Aeq,beq,lb,ub,nonlcon,options);
     [diff,SSE] = FminconFun(opt_param);
     AllParam = [diff, SSE(1), SSE(2), SSE(3), SSE(4), SSE(5), SSE(6), SSE(7), SSE(8), SSE(9) opt_param];
     dlmwrite('PosteriorA0.csv',AllParam,'-append');
     
     tsp = [0,600];
     [tp,yp]=ode15s(@eq,tsp,u0,[],opt_param);

 end

       function [diff,SSE] = FminconFun(par0)
        
        tspan=t; %we want x at every t
        [t,y]=ode15s(@eq,tspan,u0,[],par0);
        SSE(1)=0;SSE(2)=0;SSE(3)=0;SSE(4)=0;SSE(5)=0;SSE(6)=0;SSE(7)=0;
       
        if size(y(:,3),1) == noTimepoints && size(y(:,4),1) == noTimepoints && size(y(:,5),1)== noTimepoints && size(y(:,6),1)== noTimepoints && size(y(:,7),1)== noTimepoints
           diff=sum(((y(:,1) - dataPI)/mean(dataPI)).^2);
           % disp(printdiff);
        else 
           diff=9999;
        end
       end

function dz=eq(t,z,par0)
         
         r1=par0(1); 
         r1m=par0(2); 
         r2=par0(3); 
         r2m=par0(4);
         r3=par0(5); 
         r4=par0(6);
         t1=par0(7);
         t2=par0(8);
         t3=par0(9);
         t3m=par0(10);
         t4=par0(11);
         t4m=par0(12);
         t5=par0(13);
         t5m=par0(14);
         s1=par0(15);
         s2=par0(16);
         s2m=par0(17);
       
         st=(0.001*t*exp(-0.0002*t^2) + 1*tanh(0.02*t)); 
        
         dz=zeros(9,1);
     
         dz(1) = t2*z(9) - r1*z(1) + r1m*z(2) - t3*z(1) + t3m*z(6); %PI
         dz(2) = r1*z(1) - r1m*z(2) - r2*z(2) + r2m*z(3); %PI4
         dz(3) = r2*z(2) - r2m*z(3) - s2*st*z(3) + s2m*z(4) - s1*st*z(3) + t5*z(6) - t5m*z(3); %PIP2 
         dz(4) = s2*st*z(3) - s2m*z(4) - r4*z(4); %PIP3
         dz(5) = t4*z(6) - t4m*z(5) + r4*z(4); %PI34 
         dz(6) = t3*z(1) - t3m*z(6) - t5*z(6) + t5m*z(3) - t4*z(6) + t4m*z(5); %Pp
         dz(7) = s1*st*z(3) - t1*z(7); %IP3
         dz(8) = t1*z(7) - r3*z(8); %IP1
         dz(9) = r3*z(8) - t2*z(9); %Ip
    
         % t
         
         end
end
