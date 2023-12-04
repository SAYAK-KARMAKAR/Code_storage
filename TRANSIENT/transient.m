clc;
clear all;
L=[550 450];
D=[0.75 0.6];
a=[1100 900];
f=[0.01 0.012];
g=9.81;
Qo=1; % steady state discharge
Num_pipe=2; % number of pipes
Num_reach=2; %number of reach at which each pipe is divided into
Num_node=Num_pipe*(Num_reach+1);
H_ds=60; % reservior head, m
T_last=10; % time upto which computations should be done in sec
delta_t=L(1)/(a(1)*Num_reach); % in sec
Num_timestep=T_last/delta_t;
Num_timenode=Num_timestep+1;

%%
% pipe cs area
CS_area=zeros(1,2);
CS_area(1)=(pi/4)*D(1)^2;
CS_area(2)=(pi/4)*D(2)^2;

% x coordinate vector
x_pipe1=linspace(0,L(1),3);
x_pipe2=linspace(L(1),L(1)+L(2),3);
x_vector=[x_pipe1 x_pipe2];

%pipe coff Ca
Ca=zeros(1,2);
Ca(1)=(g*CS_area(1))/a(1);
Ca(2)=(g*CS_area(2))/a(2);

% coff R
R=zeros(1,2);
R(1)=f(1)/(2*D(1)*CS_area(1));
R(2)=f(2)/(2*D(2)*CS_area(2));

% Valve opening condition
Tau=zeros(1,Num_timenode);
Tau(1)=1;Tau(5)=0.9;Tau(9)=0.7;Tau(13)=0.5;Tau(17)=0.3;
Tau(21)=0.1;Tau(25)=0;
for i=0:1:5
    Tau(4*i+2)=Tau(4*i+1)+(Tau(4*i+5)-Tau(4*i+1))*(1/4);
 Tau(4*i+3)=Tau(4*i+1)+(Tau(4*i+5)-Tau(4*i+1))*(2/4);
  Tau(4*i+4)=Tau(4*i+1)+(Tau(4*i+5)-Tau(4*i+1))*(3/4);
end
%% Calculation for steady state conditions
H=zeros(Num_timenode,Num_node);
Q=zeros(Num_timenode,Num_node);
Q(1,:)=Qo;
H(1,6)=H_ds;
H(1,5)=H_ds+(f(2)*(x_vector(6)-x_vector(5))*Qo^2)/(2*D(2)*g*CS_area(2)^2);
H(1,4)=H(1,5)+(f(2)*(x_vector(5)-x_vector(4))*Qo^2)/(2*D(2)*g*CS_area(2)^2);
H(1,3)=H(1,4);
H(1,2)=H(1,3)+(f(1)*(x_vector(3)-x_vector(2))*Qo^2)/(2*D(1)*g*CS_area(1)^2);
H(1,1)=H(1,2)+(f(1)*(x_vector(2)-x_vector(1))*Qo^2)/(2*D(1)*g*CS_area(1)^2);
%junction loss and other minor losses are neglected

%% Boundary conditions
Cn=zeros(Num_timenode,Num_node);
Cp=zeros(Num_timenode,Num_node);
Cv=zeros(Num_timenode,Num_node);
while i<Num_timenode
    i=i+1;
for i=2:Num_timenode

% FOR UPSTREAM RESERVOIR***
    Cn(i,1)=Q(i-1,2)-(g*CS_area(1))/a(1)*H(i-1,2)-R(1)*delta_t*Q(i-1,2)*abs(Q(i-1,2));
    Q(i,1)=Ca(1)*H(1,1)+Cn(1);
    H(i,1)=H(1,1);

% Interior nodes
    Cp(i,2)=Q(i-1,1)+((g*CS_area(1))/a(1))*H(i-1,1)-R(1)*delta_t*Q(i-1,1)*abs(Q(i-1,1));
    Cn(i,2)=Q(i-1,3)-((g*CS_area(1))/a(1))*H(i-1,3)-R(1)*delta_t*Q(i-1,3)*abs(Q(i-1,3));
     Q(i,2)=0.5*(Cp(i,2)+Cn(i,2));
     H(i,2)=(Q(i,2)-Cn(i,2))/Ca(1);

    Cp(i,5)=Q(i-1,4)+((g*CS_area(2))/a(2))*H(i-1,4)-R(2)*delta_t*Q(i-1,4)*abs(Q(i-1,4));
    Cn(i,5)=Q(i-1,6)-((g*CS_area(2))/a(2))*H(i-1,6)-R(2)*delta_t*Q(i-1,6)*abs(Q(i-1,6));
    Q(i,5)=0.5*(Cp(i,5)+Cn(i,5));
    H(i,5)=(Q(i,5)-Cn(i,5))/Ca(2);

% Junction condition
    Cp(i,3)=Q(i-1,2)+((g*CS_area(1))/a(1))*H(i-1,2)-R(1)*delta_t*Q(i-1,2)*abs(Q(i-1,2));
    Cn(i,4)=Q(i-1,5)-((g*CS_area(2))/a(2))*H(i-1,5)-R(2)*delta_t*Q(i-1,5)*abs(Q(i-1,5));
    H(i,3)=(Cp(i,3)-Cn(i,4))/(Ca(1)+Ca(2));
    H(i,4)=H(i,3);
    Q(i,3)=Cp(i,3)-Ca(1)*H(i,3);
    Q(i,4)=Q(i,3);
% Valve BC
    Cp(i,6)=Q(i-1,5)+((g*CS_area(2))/a(2))*H(i-1,5)-R(2)*delta_t*Q(i-1,5)*abs(Q(i-1,5));
    Cv(i,6)=((Tau(i)*Qo)^2)/(Ca(2)*H(1,6));
    Q(i,6)=0.5*(-Cv(i,6)+sqrt(Cv(i,6)^2+4*Cp(i,6)*Cv(i,6)));
    H(i,6)=(Cp(i,6)-Q(i,6))/Ca(2);
end
end
Q
H
time=0:delta_t:T_last;
H_valve=transpose(H(:,6));
plot(time,H_valve);
title("Head values at downstream valve");
xlabel('time in s');
ylabel('Head in m');








































